"""mzML -> TSV scan<TAB>group, group = analyzer/activation from the Thermo filter string."""
import sys, re
def group(f):
    an = f.split()[0]
    an = {'FTMS': 'orbitrap', 'ITMS': 'ion_trap', 'ASTMS': 'astral'}.get(an, an)
    acts = re.findall(r'@([a-z]+)', f.split(' ms2 ', 1)[1])
    a = '+'.join(acts)
    a = {'hcd': 'hcd', 'cid': 'cid', 'etd': 'etd', 'etd+hcd': 'ethcd', 'etd+cid': 'etcid'}.get(a, a)
    return f'{an}/{a}'
src, out = sys.argv[1], sys.argv[2]
scan = None; n = 0
with open(src) as fh, open(out, 'w') as o:
    for line in fh:
        if '<spectrum ' in line:
            m = re.search(r'scan=(\d+)', line); scan = m.group(1) if m else None
        elif 'MS:1000512' in line and scan is not None:
            f = re.search(r'value="([^"]*)"', line).group(1)
            if ' ms2 ' in f:
                o.write(f'{scan}\t{group(f)}\n'); n += 1
            scan = None
print(src, n)
