#include <iostream>
#include <fstream>
#include <vector>
#include<queue>
#include<algorithm>
using namespace std;
const int metaSize = 12;
struct pnt {
    pnt(int i, int i1) : idx(i), lv(i1) {}
    int idx;
    int lv;
};

struct p {
    int k, v;
    p() {
        k = 0;
        v = 0;
    };

    p(int kk, int vv) : k(kk), v(vv) {};
};

bool compare(p p1, p p2) {
    return p1.k < p2.k;
}

struct SplitedBlocks {
    vector<int> left;
    vector<int> right;
};

class BPlusTree {
public:
    int blockSize, rootBID, depth;//META DATA
    int blockLength;
    long long totalNode;
    fstream file;

    BPlusTree(char *fileName) {
        //open file and init load meta data
        file.open(fileName, ios::binary | ios::out | ios::in);
        file.read(reinterpret_cast<char *>(&blockSize), 4);
        file.read(reinterpret_cast<char *>(&rootBID), 4);
        file.read(reinterpret_cast<char *>(&depth), 4);
        file.seekg(0, ios::end);
        blockLength = (blockSize) / 4;
        totalNode = (file.tellg() - 12LL) / blockSize;
    }

    ~BPlusTree() {
        file.close();
    }

    void setRoot(int idx);
    vector<int> blockToVector(int idx);
    long long insertBlock(vector<int> &v);
    long long writeBlock(vector<int> &v, int idx);


    vector<int> findLeafBlock(int key);
    vector<int> findLeafBlock(int key, int d, vector<int> curBlock, vector<int> &path);

    static void insertPairIntoLeaf(vector<int> &v, int key, int value);
    static void insertPairIntoNonLeaf(vector<int> &v, int key, int pointer);
    static void insertPairIntoBlock(vector<int> &v, int key, int pointer, bool isNonLeaf);

    SplitedBlocks splitBlock(int n, vector<int> &v, bool isLeaf) const;
    void rebalancing(vector<int> &path, vector<int> &block);
    void rebalanceNonLeafBlock(vector<int> &path, vector<int> &block);

    void insertNode(int key, int value);
    int pointSearch(int key);
    vector<p> rangeSearch(int startK, int endK);
    void printBPlusTree(ofstream& outF);
    void printNode(int idx, bool isLeaf,ofstream& outF);

};
p parseInput2p(const char* inputStr){
    int i = 0, k = 0, v = 0;
    while(inputStr[i]!=','){
        k = k * 10 + ( inputStr[i++] - '0');
    }
    i++;
    while(inputStr[i]!='\0'){
        if ('0' <= inputStr[i] && inputStr[i] <= '9') {
            v = v * 10 + ( inputStr[i] - '0');
        }
        i++;
    }
    return p(k,v);
}
int parseInput2i(const char* inputStr){
    int i = 0, k = 0;
    while(inputStr[i]!='\0'){
        if ('0' <= inputStr[i] && inputStr[i] <= '9')
            k = k * 10 + ( inputStr[i] - '0');
        i++;
    }
    return k;
}
int main(int argc, char* argv[]) {

    ofstream outF;
    ifstream inF;
    char command = argv[1][0];
    char* bin = argv[2];
    char* input;
    int blockSize=0;
    int dummy=0;
    BPlusTree *mytree;
    if (command != 'c') {
        mytree = new BPlusTree(bin);
    }
    char  inputStr[25];

    switch (command) {
        case 'c':
            blockSize = atoi(argv[3]);
            outF.open(bin, ios::binary);
            outF.write(reinterpret_cast<char*>(&blockSize), sizeof(blockSize));
            outF.write(reinterpret_cast<char*>(&dummy), sizeof(dummy));
            outF.write(reinterpret_cast<char*>(&dummy), sizeof(dummy));
            outF.close();

            break;
        case 'i':
            input = argv[3];
            inF.open(input);
            while (inF.getline( inputStr, sizeof( inputStr))) {
                p inputPair = parseInput2p(inputStr);
                mytree->insertNode(inputPair.k, inputPair.v);
            }
            inF.close();
            break;
        case 's':
            inF.open(argv[3]);
            outF.open(argv[4]);
            while (inF.getline( inputStr, sizeof( inputStr))) {
                int k = parseInput2i(inputStr);
                int v = mytree->pointSearch(k);
                outF << k << "," << v << "\n";
            }
            inF.close();
            outF.close();
            break;
        case 'r':
            inF.open(argv[3]);
            outF.open(argv[4]);
            while (inF.getline( inputStr, sizeof( inputStr))) {
                p inputPair = parseInput2p(inputStr);
                vector<p> v = mytree->rangeSearch(inputPair.k, inputPair.v);
                for (int i = 0; i < v.size(); i++) {
                    outF << v[i].k << "," << v[i].v << " / ";
                }
                outF << "\n";
            }
            inF.close();
            outF.close();
            break;
        case 'p':
            outF.open(argv[3]);
            mytree->printBPlusTree(outF);
            outF.close();
            break;
    }

    return 0;
}

void BPlusTree::setRoot(int idx) {
    if (rootBID != 0) {
        depth++;
    }
    rootBID = idx;
    file.seekg(4, ios::beg);
    file.write(reinterpret_cast<char *>(&rootBID), 4);
    file.write(reinterpret_cast<char *>(&depth), 4);
}

vector<int> BPlusTree::blockToVector(int idx) {
    vector<int> ret(blockLength);
    file.seekg(metaSize + (idx - 1) * blockSize);
    file.read(reinterpret_cast<char *>(&ret[0]), blockSize);
    return ret;
}

long long BPlusTree::writeBlock(vector<int> &v, int idx) {
    if (idx != -1) {
        file.seekp(metaSize + (idx - 1) * blockSize);
    } else {
        file.seekp(0, ios::end);
    }
    while (v.size() < blockLength) {
        v.push_back(0);
    }
    file.write(reinterpret_cast<const char *>(&v[0]), blockSize);
    if (file.tellg() != -1) {
        if (idx == -1) return ++totalNode;
        else return idx;
    }
    return -1;
}

long long BPlusTree::insertBlock(vector<int> &v) {
    return writeBlock(v, -1);
}

vector<int> BPlusTree::findLeafBlock(int key) {
    vector<int> path;
    path.push_back(rootBID);
    return findLeafBlock(key, 0, (rootBID ? blockToVector(rootBID) : vector<int>()), path);
}

vector<int> BPlusTree::findLeafBlock(int key, int d, vector<int> curBlock, vector<int> &path) {
    if (totalNode == 0)
        return vector<int>();
    if (totalNode == 1)
        return blockToVector(1);
    if (d == depth) {
        return curBlock;
    } else {
        vector<int> ret;
        for (int i = 0; i < curBlock.size() / 2; ++i) {
            int search_key = curBlock[(i << 1) + 1];
            int bId = curBlock[i << 1];
            if (key < search_key) {
                path.push_back(bId);
                return findLeafBlock(key, d + 1, blockToVector(bId), path);
            } else if (search_key == 0) {
                path.push_back(curBlock[i << 1]);
                return findLeafBlock(key, d + 1, blockToVector(curBlock[i << 1]), path);
            }
        }
        path.push_back(curBlock.back());
        return findLeafBlock(key, d + 1, blockToVector(curBlock.back()), path);
    }
}

void BPlusTree::insertPairIntoNonLeaf(vector<int> &v, int key, int pointer) {
    insertPairIntoBlock(v, key, pointer, true);
}

void BPlusTree::insertPairIntoLeaf(vector<int> &v, int key, int value) {
    insertPairIntoBlock(v, key, value, false);
}

void BPlusTree::insertNode(int key, int value) {
    vector<int> path;
    path.push_back(rootBID);
    vector<int> block = findLeafBlock(key, 0, (rootBID ? blockToVector(rootBID) : vector<int>()), path);
    insertPairIntoLeaf(block, key, value);
    if (rootBID == 0) {
        setRoot(1);
        path[0] = rootBID;
        insertBlock(block);
    }
    if (block.size() > blockLength) {
        rebalancing(path, block);
    } else if (block.size() <= blockLength) {
        int temp = block.back();
        block.pop_back();
        while (block.size() < blockLength) {
            block.push_back(0);
        }
        block.back() = (temp);

        writeBlock(block, path.back());

    }
}

void BPlusTree::rebalancing(vector<int> &path, vector<int> &block) {
    int leftLen = (blockLength + 1) / 2 + 1;
    SplitedBlocks sp = splitBlock(leftLen, block, true);
    long long idx = insertBlock(sp.right);
    sp.left.back() = idx;
    writeBlock(sp.left, path.back());
    int upperKey = sp.right[0];
    int upperPtr = idx;
    if (path.back() == rootBID) {
        vector<int> newRoot(3);
        newRoot[0] = path.back();
        newRoot[1] = upperKey;
        newRoot[2] = upperPtr;
        int rID = insertBlock(newRoot);
        setRoot(rID);
        return;
    }
    path.pop_back();
    vector<int> parentBlock = blockToVector(path.back());
    insertPairIntoNonLeaf(parentBlock, upperKey, upperPtr);
    if (parentBlock.size() > blockLength) {
        rebalanceNonLeafBlock(path, parentBlock);
    } else {
        while (parentBlock.size() < blockLength) {
            parentBlock.push_back(0);
        }
        writeBlock(parentBlock, path.back());
    }
}

void BPlusTree::rebalanceNonLeafBlock(vector<int> &path, vector<int> &block) {
    vector<int> &curBlock = block;
    while (curBlock.size() > blockLength) {
        int leftLength = (blockLength + 1) / 2;
        SplitedBlocks sp = splitBlock(leftLength, curBlock, false);
        int leftBlockIdx = path.back();
        int rightFirstKey = sp.right.front();
        writeBlock(sp.left, leftBlockIdx);
        sp.right.erase(sp.right.begin());
        int rightBlockIdx = insertBlock(sp.right);
        if (path.back() == rootBID) {
            vector<int> newRoot(blockLength);
            newRoot[0] = leftBlockIdx;
            newRoot[1] = rightFirstKey;
            newRoot[2] = rightBlockIdx;
            int rId = insertBlock(newRoot);
            setRoot(rId);
            return;
        } else {
            path.pop_back();
            int parentBlockIdx = path.back();
            vector<int> parentBlock = blockToVector(parentBlockIdx);
            insertPairIntoNonLeaf(parentBlock, rightFirstKey, rightBlockIdx);
            curBlock = parentBlock;
        }
    }
    writeBlock(curBlock, path.back());
}

void BPlusTree::insertPairIntoBlock(vector<int> &v, int key, int pointer, bool isNonLeaf) {
    vector<p> temp;
    for (int i = 0; i < v.size() / 2; ++i) {
        int keyIdx = (i << 1) + isNonLeaf;
        int valueIdx = (i << 1) + isNonLeaf + 1;
        if (v[keyIdx])
            temp.push_back(p(v[keyIdx], v[valueIdx]));
    }
    temp.push_back(p(key, pointer));
    int tempValue = -1;
    if (!isNonLeaf && v.size()) {
        tempValue = v.back();
    }
    v.resize(max(v.size(), temp.size() * 2 + 1));
    sort(temp.begin(), temp.end(), compare);
    for (int i = 0; i < temp.size(); ++i) {
        int keyIdx = (i << 1) + isNonLeaf;
        int valueIdx = (i << 1) + isNonLeaf + 1;
        v[keyIdx] = temp[i].k;
        v[valueIdx] = temp[i].v;
    }

    if (!isNonLeaf && tempValue != -1)
        v[v.size() - 1] = tempValue;
}

SplitedBlocks BPlusTree::splitBlock(int n, vector<int> &v, bool isLeaf) const {
    SplitedBlocks ret;
    for (int i = 0; i < n; ++i) {
        ret.left.push_back(v[i]);
    }
    while (ret.left.size() < blockLength) {
        ret.left.push_back(0);
    }
    for (int i = n; i < v.size() - 1 + (!isLeaf); ++i) {
        ret.right.push_back(v[i]);
    }
    while (ret.right.size() < blockLength) {
        ret.right.push_back(0);
    }
    if (isLeaf) {
        ret.right.back() = v.back();
    }
    return ret;
}

int BPlusTree::pointSearch(int key) {
    vector<int> block = findLeafBlock(key);
    for (int i = 0; i < block.size() / 2; ++i) {
        if (block[i << 1] == key) {
            return block[1 + (i << 1)];
        }
    }
    return 0;
}

vector<p> BPlusTree::rangeSearch(int startK, int endK) {
    vector<int> block = findLeafBlock(startK);
    vector<p> ret;
    bool flg = true;
    while (flg) {
        for (int i = 0; i < block.size() / 2; ++i) {
            if (block[i << 1] < startK) continue;
            if (block[i << 1] > endK) {
                flg = false;
                break;
            }
            ret.push_back(p(block[i << 1] , block[1 + (i << 1)]));
        }
        if (!block.back()) {
            break;
        } else {
            block = blockToVector(block.back());
        }
    }
    return ret;
}

void BPlusTree::printNode(int idx, bool isLeaf,ofstream& outF) {
    vector<int> block = blockToVector(idx);
    for (int i = 0; i < block.size() / 2; ++i) {
        if (block[(i << 1) + !isLeaf])
            outF << block[(i << 1) + !isLeaf] << ",";
    }
}

void BPlusTree::printBPlusTree(ofstream& outF) {
    queue<pnt> q;
    q.push(pnt(rootBID, 0));
    int level = -1;
    while (!q.empty()) {
        int cur = q.front().idx;
        int lv = q.front().lv;
        q.pop();
        if (level < lv) {
            level = lv;
            if(level){
                outF<<"\n";
            }
            outF << "[level" << level << "]\n";
        }
        printNode(cur, lv == depth,outF);
        vector<int> curBlock = blockToVector(cur);
        if (lv == depth||lv==1) continue;
        for (int i = 0; i < curBlock.size() / 2; ++i) {
            int curI = (i << 1);
            if (curBlock[curI])
                q.push(pnt(curBlock[curI], lv + 1));
        }
        if (curBlock.back())
            q.push(pnt(curBlock.back(), lv + 1));
    }
}