#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

struct DistanceMatrix {
    double** distances;
    int num_taxa;
};

struct TreeNode {
    int parent;
    int left_child;
    int right_child;
    int parent_dist;
    int left_dist;
    int right_dist;
    int num_descendants;
    int label;
};

struct Tree {
    struct TreeNode* TreeNodeInstance; // = (struct TreeNode*) malloc(sizeof(struct TreeNode));
    int root;
};

struct Tree* create_Tree(int num_taxa) {
    struct Tree* TreeInstance = (struct Tree*) malloc(sizeof(struct Tree));
    TreeInstance->TreeNodeInstance = (struct TreeNode*) malloc(sizeof(struct TreeNode*) * (num_taxa+1));
    TreeInstance->root = num_taxa;

    TreeInstance->TreeNodeInstance[TreeInstance->root].parent = -1;
    TreeInstance->TreeNodeInstance[TreeInstance->root].left_child = -1;
    TreeInstance->TreeNodeInstance[TreeInstance->root].right_child = -1;
    TreeInstance->TreeNodeInstance[TreeInstance->root].num_descendants = 0;

    for (int i=0; i<num_taxa; i++) {
        TreeInstance->TreeNodeInstance[i].parent = -1;
        TreeInstance->TreeNodeInstance[i].left_child = -1;
        TreeInstance->TreeNodeInstance[i].right_child = -1;
        TreeInstance->TreeNodeInstance[i].num_descendants = 0;
    }
}

struct DistanceMatrix* create_distance_matrix(int num_taxa) {
    struct DistanceMatrix* dm = (struct DistanceMatrix*) malloc(sizeof(struct DistanceMatrix));
    dm->num_taxa = num_taxa;
    dm->distances = (double**) malloc(sizeof(double*) * num_taxa);
    for (int i = 0; i < num_taxa; i++) {
        dm->distances[i] = (double*) malloc(sizeof(double) * num_taxa);
        memset(dm->distances[i], 0, sizeof(double) * num_taxa);
    }
    return dm;
}

void destroy_distance_matrix(struct DistanceMatrix* dm) {
    for (int i = 0; i < dm->num_taxa; i++) {
        free(dm->distances[i]);
    }
    free(dm->distances);
    free(dm);
}

/// @brief Calculates distance (ie total nucleotide difference/ nucleotides compared) between 2 sequences
/// @param seq1 
/// @param seq2 
/// @return distance (as described in biref)
double distance(char* seq1, char* seq2) {
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);
    double diffs = 0;
    int len = fmin(len1, len2);
    for (int i = 0; i < len; i++) {
        if (seq1[i] != seq2[i]) {
            diffs++;
        }
    }
    return diffs / len;
}

/// @brief Add distances to Distance Matrix: For an array of genetic sequences, calculate distances for all combinations and store them in dm (in both indexes reflecting across diagonal)
/// @param dm 
/// @param seqs 
/// @param num_taxa 
void compute_distances(struct DistanceMatrix* dm, char** seqs, int num_taxa) {
    for (int i = 0; i < num_taxa; i++) {
        for (int j = i + 1; j < num_taxa; j++) {
            dm->distances[i][j] = dm->distances[j][i] = distance(seqs[i], seqs[j]);
        }
    }
}


/// @brief Calculate D'ij = (n-2)*Dij - TotalDistance(i) - TotalDistance(j)
/// @param dm 
/// @param i 
/// @param j 
/// @return 
double compute_q(struct DistanceMatrix* dm, int i, int j) {
    int num_taxa = dm->num_taxa;
    double sum_i = 0;
    double sum_j = 0;
    for (int k = 0; k < num_taxa; k++) {
        sum_i += dm->distances[i][k];
        sum_j += dm->distances[j][k];
    }
    return (num_taxa - 2) * dm->distances[i][j] - sum_i - sum_j;
}


/// @brief Calculate indexes with minimum D' and store in array variable pair
/// @param dm 
/// @param pair 
/// @param num_taxa 
/// @return 
double find_closest_pair(struct DistanceMatrix* dm, int* pair, int num_taxa) {
    double min_distance = INT_MAX;
    for (int i = 0; i < num_taxa; i++) {
        for (int j = i + 1; j < num_taxa; j++) {
            double q = compute_q(dm, i, j);
            if (q < min_distance) {
                min_distance = q;
                pair[0] = i;
                pair[1] = j;
            }
        }
    }
    return min_distance / (num_taxa - 2);
}

/// @brief Takes the 2 nodes identied and combine them to create new node and add to the Tree
int add_node(struct TreeNode* tree, int parent, int left_child, int right_child) {
    int new_node;
    for (int i = 0; i < tree->num_descendants; i++) {
        if (tree[i].parent == -1) {
            new_node = i;
            break;
        }
    }

    tree[new_node].parent = parent;
    tree[new_node].left_child = left_child;
    tree[new_node].right_child = right_child;
    tree[new_node].num_descendants = tree[left_child].num_descendants + tree[right_child].num_descendants;

    return new_node;
}

/// @brief 
// Not sure if this will be used though
/// @param tree 
/// @param dm 
/// @param node 
/// @param distance 
void remove_node(struct TreeNode* tree, int node) {
    tree[node].parent = -1;
    tree[node].left_child = -1;
    tree[node].right_child = -1;
    tree[node].num_descendants = -1;
}

/// Update the distances in the matrix
void update_distances(struct DistanceMatrix* dm, int i, int j, double* dists) {
    int num_taxa = dm->num_taxa;
    for (int k = 0; k < num_taxa; k++) {
        if (k != i && k != j) {
            double d_ik = dm->distances[i][k];
            double d_jk = dm->distances[j][k];
            double new_distance = (d_ik + d_jk - dists[i] - dists[j]) / 2;
                dm->distances[i][k] = dm->distances[k][i] = new_distance;
                dm->distances[j][k] = dm->distances[k][j] = new_distance;
        }
    }
}

/// Not sure if we need this right now
void compute_distances_from_root(struct TreeNode* tree, struct DistanceMatrix* dm, int node, double distance) {
    if (tree[node].parent == -1) {
        return;
    }
        int parent = tree[node].parent;
        dm->distances[node][parent] = dm->distances[parent][node] = distance;
        compute_distances_from_root(tree, dm, parent, distance + dm->distances[node][parent]);
}

void nj(struct DistanceMatrix* dm, struct Tree* tree, int num_taxa) {
    if (num_taxa == 2) {
        tree->TreeNodeInstance[0].parent = tree->root;
        tree->TreeNodeInstance[0].left_child = 1;
        tree->TreeNodeInstance[0].right_child = -1;
        tree->TreeNodeInstance[0].num_descendants = 2;
        tree->TreeNodeInstance[1].parent = 0;
        tree->TreeNodeInstance[1].left_child = -1;
        tree->TreeNodeInstance[1].right_child = -1;
        tree->TreeNodeInstance[1].num_descendants = 1;
        return;
    }
    
    int* closest_pair = (int*) malloc(sizeof(int) * 2);
    double* dists = (double*) malloc(sizeof(double) * num_taxa);
    memset(dists, 0, sizeof(double) * num_taxa);
    find_closest_pair(dm, closest_pair, num_taxa);
    
    int i = closest_pair[0];
    int j = closest_pair[1];
    
    double q = compute_q(dm, i, j);
    double delta_ij = (dists[i] - dists[j] + q) / 2;
    
    int new_node = add_node(tree, -1, i, j);
    remove_node(tree, i);
    remove_node(tree, j);
    
    update_distances(dm, i, j, dists);
    compute_distances_from_root(tree, dm, new_node, delta_ij);
    nj(dm, tree, num_taxa - 1);
    
    free(closest_pair);
    free(dists);
}





int main() {
    
    int num_taxa = 4;
    char** seqs = (char**) malloc(sizeof(char*) * num_taxa);
    
    
    for (int i = 0; i < num_taxa; i++) 
    {
        seqs[i] = (char*) malloc(sizeof(char) * 5);
    }

    strcpy(seqs[0], "ACTGATCGTAGCTAGCTAGCCTCGCTAGCT");
    strcpy(seqs[1], "AGTGGATCGTAGCTAGCCCCCGTAGCTGGC");
    strcpy(seqs[2], "GCTGACTCTGCTCGTATCGATCGTACGTTG");
    strcpy(seqs[3], "AGCCGTAGCTCGTTGTATCGATGCTTTGAG");
    
    struct DistanceMatrix* dm = create_distance_matrix( num_taxa);
    compute_distances(dm, seqs, num_taxa);

    printf("Print the Distance Matrix\n");
    for (int i=0 ; i<num_taxa; i++){
        for (int j=0; j<num_taxa; j++) {
            printf("%d ", dm->distances[i][j]);
        }
        printf("\n");
    }
    //struct Tree* TreeInstance = create_Tree(num_taxa);
    //nj(dm, TreeInstance, numa_taxa);

return 0;
}