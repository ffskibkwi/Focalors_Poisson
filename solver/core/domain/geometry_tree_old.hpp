#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cstdint>
#include <iostream>
#include "core/base/location_boundary.h"
#include "core/domain/domain2d.h"

// Currently the main project uses 2D, but the tree structure and building algorithm are independent of dimension.
//DomainT = Domain2DUniform / (Domain2DNonuniform) / Domain3DUniform / (Domain3DNonuniform)

template <typename DomainT>
struct GeometryTreeNode {
    DomainT* domain;
    std::vector<GeometryTreeNode*> children;
    explicit GeometryTreeNode(DomainT* d) : domain(d) {}
};

template <typename DomainT>
class TreeBuilder {
public:
    using Adjacency = std::unordered_map<DomainT*, std::unordered_map<LocationType, DomainT*>>;

    GeometryTreeNode<DomainT>* buildOptimalTree(const Adjacency& adjacency)
    {
        if (adjacency.empty()) {
            return nullptr;
        }

        // Collect all nodes and calculate the degree
        std::unordered_map<DomainT*, int> degrees;
        std::unordered_set<DomainT*> nodes;
        for (const auto& pair : adjacency) {
            nodes.insert(pair.first);
            for (const auto& nPair : pair.second) {
                nodes.insert(nPair.second);
            }
        }
        for (DomainT* node : nodes) {
            if (adjacency.count(node)) {
                degrees[node] = static_cast<int>(adjacency.at(node).size());
            } else {
                degrees[node] = 0;
            }
        }

        // Find the graph center
        std::queue<DomainT*> leaves;
        for (DomainT* node : nodes) {
            if (degrees[node] <= 1) {
                leaves.push(node);
            }
        }
        int remainingNodes = static_cast<int>(nodes.size());
        while (remainingNodes > 2) {
            int leavesCount = static_cast<int>(leaves.size());
            if (leavesCount == 0) break; // If it is not a tree structure, prevent infinite loop
            remainingNodes -= leavesCount;
            for (int i = 0; i < leavesCount; ++i) {
                DomainT* leaf = leaves.front();
                leaves.pop();
                if (adjacency.count(leaf)) {
                    for (const auto& p : adjacency.at(leaf)) {
                        DomainT* neighbor = p.second;
                        degrees[neighbor]--;
                        if (degrees[neighbor] == 1) {
                            leaves.push(neighbor);
                        }
                    }
                }
            }
        }

        // Choose the unique root
        if (leaves.empty()) {
            if (!nodes.empty()) {
                leaves.push(*nodes.begin());
            } else {
                return nullptr;
            }
        }
        DomainT* rootDomain = leaves.front();
        leaves.pop();
        if (!leaves.empty()) {
            DomainT* secondCenter = leaves.front();
            if (reinterpret_cast<uintptr_t>(secondCenter) < reinterpret_cast<uintptr_t>(rootDomain)) {
                rootDomain = secondCenter;
            }
        }

        // Recursively build the tree
        return buildSubtree(rootDomain, nullptr, adjacency);
    }

private:
    GeometryTreeNode<DomainT>* buildSubtree(
        DomainT* currentDomain,
        DomainT* parentDomain,
        const Adjacency& adjacency)
    {
        // Use an explicit visited set to avoid infinite recursion caused by loops in general graphs.
        std::unordered_set<DomainT*> visited;
        return buildSubtreeImpl(currentDomain, parentDomain, adjacency, visited);
    }

    GeometryTreeNode<DomainT>* buildSubtreeImpl(
        DomainT* currentDomain,
        DomainT* parentDomain,
        const Adjacency& adjacency,
        std::unordered_set<DomainT*>& visited)
    {
        if (visited.count(currentDomain)) return nullptr;
        visited.insert(currentDomain);

        GeometryTreeNode<DomainT>* node = new GeometryTreeNode<DomainT>(currentDomain);
        if (adjacency.count(currentDomain)) {
            for (const auto& pair : adjacency.at(currentDomain)) {
                DomainT* neighbor = pair.second;
                if (neighbor == parentDomain) continue;
                if (auto* child = buildSubtreeImpl(neighbor, currentDomain, adjacency, visited)) {
                    node->children.push_back(child);
                }
            }
        }
        return node;
    }
};

template <typename DomainT>
void deleteTree(GeometryTreeNode<DomainT>* root)
{
    if (!root) return;
    for (auto* child : root->children) {
        deleteTree(child);
    }
    delete root;
}

template <typename DomainT>
void printTree(GeometryTreeNode<DomainT>* node, int depth = 0)
{
    if (!node) return;
    for (int i = 0; i < depth; ++i) std::cout << "  ";
    std::cout << "- " << node->domain->name << std::endl;
    for (auto* child : node->children) {
        printTree(child, depth + 1);
    }
}

// Only for 2D
using GeometryTreeNode2D = GeometryTreeNode<Domain2DUniform>;
using TreeBuilder2D = TreeBuilder<Domain2DUniform>;