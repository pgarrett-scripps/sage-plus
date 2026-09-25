use crate::cloud_store::{CloudError, CloudProvider};

impl CloudProvider {
    /// Recognize cloud schemes even without the `cloud` feature, so cloud URIs
    /// fail with `FeatureNotEnabled` instead of being treated as local paths.
    pub(crate) fn parse(url: impl AsRef<str>) -> Option<Self> {
        let url = url.as_ref();
        let scheme = url.find("://").map_or("", |pos| &url[..pos]);
        match scheme {
            "s3" | "s3a" => Some(Self::S3),
            "az" | "adl" | "azure" | "abfs" | "abfss" => Some(Self::Azure),
            "gs" => Some(Self::Gcs),
            _ => None,
        }
    }
}

#[derive(Debug)]
pub(crate) struct CloudObject {}

impl CloudObject {
    pub(crate) fn new(_url: impl AsRef<str>) -> Result<Self, CloudError> {
        Err(CloudError::FeatureNotEnabled)
    }

    #[allow(clippy::len_without_is_empty)]
    pub(crate) fn len(&self) -> Result<usize, CloudError> {
        Err(CloudError::FeatureNotEnabled)
    }

    pub(crate) fn children(&self) -> Result<Vec<String>, CloudError> {
        Err(CloudError::FeatureNotEnabled)
    }

    pub(crate) fn range(
        &self,
        _range: impl std::ops::RangeBounds<usize>,
    ) -> Result<Vec<u8>, CloudError> {
        Err(CloudError::FeatureNotEnabled)
    }

    pub(crate) fn download_to(
        &self,
        _dest: &std::path::Path,
    ) -> Result<(), CloudError> {
        Err(CloudError::FeatureNotEnabled)
    }

    pub(crate) fn upload_bytes(
        &self,
        _bytes: Vec<u8>,
    ) -> Result<(), CloudError> {
        Err(CloudError::FeatureNotEnabled)
    }

    pub(crate) fn upload_from(
        &self,
        _src: impl AsRef<std::path::Path>,
    ) -> Result<(), CloudError> {
        Err(CloudError::FeatureNotEnabled)
    }
}
