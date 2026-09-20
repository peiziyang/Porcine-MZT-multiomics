"""Check mofapy2 internals for factor extraction."""
import mofapy2, os
m_dir = os.path.dirname(mofapy2.__file__)
# Check entry_point for node access
ep_file = os.path.join(m_dir, "run", "entry_point.py")
if os.path.exists(ep_file):
    with open(ep_file) as f:
        for line in f:
            if "def get_" in line or "def getZ" in line or "def getW" in line or ".node" in line:
                print(line.strip())
# Check node base class
node_file = os.path.join(m_dir, "nodes", "node.py")
if os.path.exists(node_file):
    with open(node_file) as f:
        for line in f:
            if "def getExpect" in line or "def getFactors" in line or "def getWeights" in line:
                print(line.strip())
print("\n--- nodes directory ---")
for f in sorted(os.listdir(os.path.join(m_dir, "nodes"))):
    if f.endswith(".py"):
        print(f)
