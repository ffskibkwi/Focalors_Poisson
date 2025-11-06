#pragma once
#include "base/domain/domain2d.h"
#include "base/location_boundary.h"
#include <cstdint>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace TreeUtils
{
    template<typename DomainT>
    using Adjacency = std::unordered_map<DomainT*, std::unordered_map<LocationType, DomainT*>>;
    /**
     * @brief 从图的邻接表中找到最优的根节点（图的中心）。
     * @tparam DomainT 节点的类型。
     * @param adjacency 图的邻接表。
     * @return 指向最优根节点的指针；如果图为空，则返回 nullptr。
     */
    template<typename DomainT>
    DomainT* findOptimalRoot(const Adjacency<DomainT>& adjacency)
    {
        if (adjacency.empty())
        {
            return nullptr;
        }

        std::unordered_set<DomainT*>      nodes;
        std::unordered_map<DomainT*, int> degrees;
        std::queue<DomainT*>              leaves;

        // Init nodes
        for (const auto& pair : adjacency)
            nodes.insert(pair.first);

        // Init degrees
        for (DomainT* node : nodes)
            degrees[node] = static_cast<int>(adjacency.at(node).size());

        // Init leaves
        for (DomainT* node : nodes)
            if (degrees[node] <= 1)
                leaves.push(node);

        // Strip leaves layer by layer
        int remainingNodes = static_cast<int>(nodes.size());
        while (remainingNodes > 2)
        {
            // Validate
            int leavesCount = static_cast<int>(leaves.size());
            if (leavesCount == 0)
                break;

            // Update counter
            remainingNodes -= leavesCount;

            for (int i = 0; i < leavesCount; ++i)
            {
                // Strip
                DomainT* leaf = leaves.front();
                leaves.pop();
                for (const auto& p : adjacency.at(leaf))
                {
                    DomainT* neighbor = p.second;
                    if (--degrees[neighbor] == 1)
                    {
                        // Striping will introduce new leaves
                        leaves.push(neighbor);
                    }
                }
            }
        }

        // If graph has two nodes and has no leaves,
        // it is a cycle
        if (leaves.empty())
            return nodes.empty() ? nullptr : *nodes.begin();

        // Choose root from remaining two or one
        DomainT* rootDomain = leaves.front();
        leaves.pop();
        // If remaining two
        if (!leaves.empty())
        {
            // Choose the smaller address
            DomainT* secondCenter = leaves.front();
            if (reinterpret_cast<uintptr_t>(secondCenter) < reinterpret_cast<uintptr_t>(rootDomain))
            {
                rootDomain = secondCenter;
            }
        }
        return rootDomain;
    }
    // 内部递归辅助函数
    template<typename DomainT>
    void buildTreeMapRecursive(DomainT*                                                                  current,
                               DomainT*                                                                  parent,
                               const Adjacency<DomainT>&                                                 adjacency,
                               std::unordered_map<DomainT*, std::unordered_map<LocationType, DomainT*>>& tree_map)
    {
        if (!current)
            return;
        tree_map.try_emplace(current); // 确保当前节点作为键存在
        if (adjacency.count(current))
        {
            for (const auto& pair : adjacency.at(current))
            {
                const LocationType childDir = pair.first;
                DomainT*           neighbor = pair.second;
                // 只需要检查父节点即可防止回溯，因为输入保证是树状结构
                if (neighbor == parent)
                {
                    continue;
                }
                tree_map[current][childDir] = neighbor;
                buildTreeMapRecursive(neighbor, current, adjacency, tree_map);
            }
        }
    }
    /**
     * @brief 从一个给定的根节点开始，构建树的哈希表表示。
     * @tparam DomainT 节点的类型。
     * @param root 开始构建树的根节点。
     * @param adjacency 图的邻接表。
     * @return 一个表示树结构的哈希表。
     */
    template<typename DomainT>
    std::unordered_map<DomainT*, std::unordered_map<LocationType, DomainT*>>
    buildTreeMapFromRoot(DomainT* root, const Adjacency<DomainT>& adjacency)
    {
        if (!root)
        {
            return {};
        }
        std::unordered_map<DomainT*, std::unordered_map<LocationType, DomainT*>> tree_map;
        buildTreeMapRecursive<DomainT>(root, nullptr, adjacency, tree_map);
        return tree_map;
    }

    /**
     * @brief 从树的哈希表中构建父节点哈希表。
     * @tparam DomainT 节点的类型。
     * @param tree_map 树的哈希表。
     * @return 一个表示父节点哈希表。
     */
    template<typename DomainT>
    std::unordered_map<DomainT*, std::pair<LocationType, DomainT*>>
    buildParentMapFromTree(std::unordered_map<DomainT*, std::unordered_map<LocationType, DomainT*>>& tree_map)
    {
        std::unordered_map<DomainT*, std::pair<LocationType, DomainT*>> parent_map;
        for (auto& [domain, children] : tree_map)
        {
            for (auto& [location, child] : children)
            {
                parent_map[child] = std::make_pair(opposite(location), domain);
            }
        }
        return parent_map;
    }

    inline const char* locationTypeToString(LocationType t)
    {
        switch (t)
        {
            case LocationType::Left:
                return "Left";
            case LocationType::Right:
                return "Right";
            case LocationType::Down:
                return "Down";
            case LocationType::Up:
                return "Up";
            case LocationType::Front:
                return "Front";
            case LocationType::Back:
                return "Back";
            default:
                return "Unknown";
        }
    }

    // --- 打印辅助函数 ---
    template<typename DomainT>
    void printTreeMap(DomainT*                                                                        currentNode,
                      const std::unordered_map<DomainT*, std::unordered_map<LocationType, DomainT*>>& tree_map,
                      int                                                                             depth = 0,
                      bool printNodeLine                                                                    = true)
    {
        if (!currentNode)
            return;
        if (printNodeLine)
        {
            for (int i = 0; i < depth; ++i)
                std::cout << "  ";
            std::cout << "- " << currentNode->name << std::endl;
        }
        if (tree_map.count(currentNode))
        {
            int dirIndent = depth + (printNodeLine ? 1 : 0);
            for (const auto& kv : tree_map.at(currentNode))
            {
                for (int i = 0; i < dirIndent; ++i)
                    std::cout << "  ";
                std::cout << TreeUtils::locationTypeToString(kv.first) << " -> " << kv.second->name << std::endl;
                printTreeMap(kv.second, tree_map, dirIndent + 1, false);
            }
        }
    }
} // namespace TreeUtils