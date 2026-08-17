use rmcp::{
    model::CallToolRequestParams,
    transport::{ConfigureCommandExt, TokioChildProcess},
    ServiceExt,
};
use tokio::process::Command;
use uuid::Uuid;

fn arguments(value: serde_json::Value) -> serde_json::Map<String, serde_json::Value> {
    serde_json::from_value(value).expect("tool arguments must be an object")
}

#[tokio::test]
async fn exposes_sage_tools_over_stdio() -> anyhow::Result<()> {
    let root = std::env::temp_dir().join(format!("sage-mcp-stdio-{}", Uuid::new_v4()));
    std::fs::create_dir_all(&root)?;
    let executable = env!("CARGO_BIN_EXE_sage-mcp");
    let transport = TokioChildProcess::new(Command::new(executable).configure(|command| {
        command.arg("--root").arg(&root);
    }))?;
    let client = ().serve(transport).await?;

    let tools = client.list_all_tools().await?;
    let names = tools
        .iter()
        .map(|tool| tool.name.as_ref())
        .collect::<Vec<_>>();
    assert!(names.contains(&"validate_config"));
    assert!(names.contains(&"inspect_config"));
    assert!(names.contains(&"estimate_search"));
    assert!(names.contains(&"start_search"));
    assert!(names.contains(&"get_job_status"));
    assert!(names.contains(&"cancel_search"));
    assert!(names.contains(&"summarize_run"));
    assert!(names.contains(&"query_results"));

    let result = client
        .call_tool(CallToolRequestParams::new("list_jobs"))
        .await?;
    assert_eq!(result.is_error, Some(false));
    assert_eq!(result.structured_content, Some(serde_json::json!([])));

    client.cancel().await?;
    std::fs::remove_dir_all(root)?;
    Ok(())
}

#[tokio::test]
async fn runs_fixture_search_through_mcp_tools() -> anyhow::Result<()> {
    let root = std::env::temp_dir().join(format!("sage-mcp-search-{}", Uuid::new_v4()));
    let tests = root.join("tests");
    std::fs::create_dir_all(&tests)?;
    let source = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests");
    for name in ["config.json", "Q99536.fasta", "LQSRPAAPPAPGPGQLTLR.mzML"] {
        std::fs::copy(source.join(name), tests.join(name))?;
    }

    let executable = env!("CARGO_BIN_EXE_sage-mcp");
    let transport = TokioChildProcess::new(Command::new(executable).configure(|command| {
        command.arg("--root").arg(&root);
    }))?;
    let client = ().serve(transport).await?;

    let validation = client
        .call_tool(
            CallToolRequestParams::new("validate_config").with_arguments(arguments(
                serde_json::json!({ "config_path": "tests/config.json" }),
            )),
        )
        .await?;
    assert_eq!(validation.is_error, Some(false));
    assert_eq!(
        validation.structured_content.as_ref().unwrap()["valid"],
        true
    );

    let inspection = client
        .call_tool(
            CallToolRequestParams::new("inspect_config").with_arguments(arguments(
                serde_json::json!({ "config_path": "tests/config.json" }),
            )),
        )
        .await?;
    assert_eq!(inspection.is_error, Some(false));
    assert_eq!(
        inspection.structured_content.as_ref().unwrap()["validation"]["valid"],
        true
    );

    let estimate = client
        .call_tool(
            CallToolRequestParams::new("estimate_search").with_arguments(arguments(
                serde_json::json!({ "config_path": "tests/config.json" }),
            )),
        )
        .await?;
    assert_eq!(estimate.is_error, Some(false));
    assert!(
        estimate.structured_content.as_ref().unwrap()["modified_peptides"]
            .as_u64()
            .unwrap()
            > 0
    );

    let started = client
        .call_tool(
            CallToolRequestParams::new("start_search").with_arguments(arguments(
                serde_json::json!({
                    "config_path": "tests/config.json",
                    "approved": true,
                    "batch_size": 1
                }),
            )),
        )
        .await?;
    assert_eq!(started.is_error, Some(false));
    let job_id = started.structured_content.as_ref().unwrap()["job_id"]
        .as_str()
        .unwrap()
        .to_owned();

    let mut final_status = None;
    for _ in 0..100 {
        let status = client
            .call_tool(
                CallToolRequestParams::new("get_job_status")
                    .with_arguments(arguments(serde_json::json!({ "job_id": job_id }))),
            )
            .await?;
        let body = status.structured_content.unwrap();
        match body["status"].as_str() {
            Some("completed" | "failed" | "cancelled") => {
                final_status = Some(body);
                break;
            }
            _ => tokio::time::sleep(std::time::Duration::from_millis(20)).await,
        }
    }
    let status = final_status.expect("fixture search did not finish in time");
    assert_eq!(status["status"], "completed", "{status:#}");
    assert!(status["summary"]["output_paths"].as_array().unwrap().len() >= 2);

    let events = client
        .call_tool(
            CallToolRequestParams::new("get_job_events").with_arguments(arguments(
                serde_json::json!({ "job_id": job_id, "after_sequence": 0 }),
            )),
        )
        .await?;
    let events = events.structured_content.unwrap();
    assert!(events
        .as_array()
        .unwrap()
        .iter()
        .any(|event| { event["event"] == "job_completed" }));

    let summary = client
        .call_tool(
            CallToolRequestParams::new("summarize_run")
                .with_arguments(arguments(serde_json::json!({ "job_id": job_id }))),
        )
        .await?;
    assert_eq!(summary.is_error, Some(false));
    assert_eq!(
        summary.structured_content.as_ref().unwrap()["status"],
        "completed"
    );

    let query = client
        .call_tool(
            CallToolRequestParams::new("query_results").with_arguments(arguments(
                serde_json::json!({
                    "job_id": job_id,
                    "dataset": "psms",
                    "limit": 5
                }),
            )),
        )
        .await?;
    assert_eq!(query.is_error, Some(false));
    assert!(
        query.structured_content.as_ref().unwrap()["returned_rows"]
            .as_u64()
            .unwrap()
            > 0
    );

    client.cancel().await?;
    std::fs::remove_dir_all(root)?;
    Ok(())
}
