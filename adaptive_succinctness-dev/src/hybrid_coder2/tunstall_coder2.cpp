#include "tunstall_coder2.hpp"

using namespace std;

namespace tunstall_coder2 {
    void traverse(vector<tree_node>& tree, uint64_t curnode, uint64_t& curindex, uint64_t sigma,
              vector<uint64_t>& currcode, Tunstall_coder_test& tunstall) {
    uint64_t i, cursum;

    for (i = 0; i < sigma; ++i, ++curnode) {
        currcode.push_back(i);

        if (tree[curnode].second == -1) {
            tree[curnode].first = curindex;  // dictionary index of this code
            cursum = 0;
            for (uint64_t j = 0; j < currcode.size(); ++j) {
                cursum += tunstall.map_table[currcode[j]];
                tunstall.D[curindex].push_back(cursum);
            }
            curindex++;
        } else {
            tree[curnode].first = curindex;
            cursum = 0;
            for (uint64_t j = 0; j < currcode.size(); ++j) {
                cursum += tunstall.map_table[currcode[j]];
                tunstall.D[curindex].push_back(cursum);
            }
            curindex++;
            traverse(tree, tree[curnode].second, curindex, sigma, currcode);  // goes to children of curnode
        }
        currcode.pop_back();
    }
}


void build_tunstall(vector<uint32_t>& seq, uint32_t block_size, uint64_t D_size_init, Tunstall_coder_test& tunstall) {
    uint64_t i;
    tunstall.D_size = D_size_init;

    tunstall.bSize = block_size;

    map<uint32_t, uint32_t> alphabet_map;

    for (i = 0; i < seq.size(); ++i)
        alphabet_map[seq[i]] = 1;

    tunstall.sigma = alphabet_map.size();

    i = 0;
    for (map<uint32_t, uint32_t>::iterator it = alphabet_map.begin(); it != alphabet_map.end(); ++it) {
        it->second = i;
        i++;
        tunstall.map_table.push_back(it->first);
    }

    vector<uint64_t> freq_table(tunstall.sigma, 0);

    for (i = 0; i < seq.size(); ++i) {
        ++freq_table[alphabet_map[seq[i]]];
    }

    uint64_t Pmin = freq_table[0];

    for (i = 1; i < freq_table.size(); ++i)
        if (freq_table[i] < Pmin)
            Pmin = freq_table[i];

    priority_queue<heap_node> H;

    vector<tree_node> tree(tunstall.sigma);

    // tree initially contains sigma elements,
    // and they have no children in the tree (indicated with -1)

    for (i = 0; i < tunstall.sigma; ++i) {
        tree[i] = tree_node(0, -1);
        H.push(heap_node((float)freq_table[i] / seq.size(), i));
    }

    // now, tree nodes are expanded according to their probabilities

    float probMinLeaf = (float)Pmin / freq_table.size();
    uint dictionary_size = tunstall.sigma;
    while (dictionary_size + tunstall.sigma <= tunstall.D_size) {
        heap_node N = H.top();
        H.pop();
        tree[N.second].second = tree.size();  // pointer to the children

        for (i = 0; i < tunstall.sigma; ++i) {
            pair<float, uint64_t> p(N.first * ((float)freq_table[i] / seq.size()), tree.size());
            H.push(p);

            tree.push_back(tree_node(0, -1));
        }
        dictionary_size = dictionary_size + tunstall.sigma;
    }

    // Now, traverse the tree to store the codes into the dictionary
    uint64_t curindex = 0;

    tunstall.D = std::vector<vector<uint64_t>>(tunstall.D_size);

    vector<uint64_t> currcode;

    traverse(tree, 0, curindex, tunstall.sigma, currcode);

    uint64_t curnode = 0;
    uint64_t prefix_sum = 0, nelems_block = 0;

    blockElement bElem;

    bElem.prefix_sum = 0;
    bElem.starting_position = 0;

    tunstall.block.push_back(bElem);

    for (i = 0; i < seq.size(); ++i) {
        prefix_sum += seq[i];
        nelems_block++;
        if (tree[curnode + alphabet_map[seq[i]]].second != -1) {
            if (nelems_block == block_size) {
                tunstall.compressed_seq.push_back(tree[curnode + alphabet_map[seq[i]]].first);
                bElem.prefix_sum = prefix_sum;
                bElem.starting_position = tunstall.compressed_seq.size();
                tunstall.block.push_back(bElem);
                nelems_block = 0;
                curnode = 0;
                nelems_block = 0;
            } else
                curnode = tree[curnode + alphabet_map[seq[i]]].second;

        } else {
            tunstall.compressed_seq.push_back(tree[curnode + alphabet_map[seq[i]]].first);
            curnode = 0;  // go back to the Tunstall tree root again
            if (nelems_block == block_size) {
                bElem.prefix_sum = prefix_sum;
                bElem.starting_position = tunstall.compressed_seq.size();
                tunstall.block.push_back(bElem);
                nelems_block = 0;
            }
        }
    }

    if (nelems_block > 0)
        tunstall.compressed_seq.push_back(tree[curnode + alphabet_map[seq[seq.size() - 1]]].first);
}

uint32_t decode(Tunstall_coder_test& tunstall, uint64_t i) {
    uint64_t b = i / tunstall.bSize;

    uint64_t sum = tunstall.block[b].prefix_sum;
    uint64_t p = tunstall.block[b].starting_position;

    uint64_t j, size, nDecode = i % tunstall.bSize + 1;

    for (j = 0; j <= nDecode; ++p) {
        size = tunstall.D[tunstall.compressed_seq[p]].size();
        if (j + size < nDecode) {
            sum += tunstall.D[tunstall.compressed_seq[p]][size - 1];
            j += size;
        } else {
            sum += tunstall.D[tunstall.compressed_seq[p]][nDecode - j - 1];
            break;
        }
    }

    return sum;
}

// Tunstall dictionary size, in bytes
uint64_t dict_size(Tunstall_coder_test& tunstall) {
    uint64_t i;
    uint64_t total_size = 0;
    for (i = 0; i < tunstall.D.size(); ++i) {
        total_size += sizeof(uint16_t) * D[i].size();
    }
    return total_size;
}

// compressed size, in bytes
uint64_t size(Tunstall_coder_test& tunstall) {
    return dict_size() + tunstall.compressed_seq.size() * sizeof(uint16_t) + tunstall.block.size() * sizeof(blockElement);
}

uint64_t nCodewords(Tunstall_coder_test& tunstall) {
    return tunstall.compressed_seq.size();
}
}