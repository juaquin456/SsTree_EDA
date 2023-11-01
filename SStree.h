#ifndef SSTREE_H
#define SSTREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
#include <limits>
#include <fstream>
//#include "SFML/Graphics.hpp"

#include "params.h"
#include "Point.h"

inline auto comp = [](std::pair<NType,std::pair<Point, std::string>> t1, std::pair<NType,std::pair<Point, std::string>> t2){return t1.first < t2.first;};
using pq_type = std::priority_queue<std::pair<NType, std::pair<Point, std::string>>, std::vector<std::pair<NType,std::pair<Point, std::string>>>, decltype(comp)>;



struct Pair {
    Point point;
    NType distance;

    Pair(const Point& p, NType d) : point(p), distance(d) {}
};
struct Comparator {
    bool operator()(const Pair& a, const Pair& b) const {
        return a.distance < b.distance; // max-heap basado en distancia
    }
};

class SsNode {
private:
    NType varianceAlongDirection(const std::vector<Point>& centroids, size_t direction) const;
    size_t minVarianceSplit(size_t coordinateIndex);

public:
    size_t D;
    SsNode(size_t dim) : D(dim) {}
//    SsNode() = default;
    virtual ~SsNode() = default;
    virtual int count() const = 0;
    Point centroid;
    NType radius;
    NType r_min;
    SsNode* parent = nullptr;
//    virtual void draw(sf::RenderWindow& window) const = 0;
    virtual bool isLeaf() const = 0;
    virtual std::vector<Point> getEntriesCentroids() const = 0;
    virtual void sortEntriesByCoordinate(size_t coordinateIndex) = 0;
    virtual std::pair<SsNode*, SsNode*> split() = 0;
    virtual bool intersectsPoint(const Point& point) const {
        return distance(this->centroid, point) <= this->radius;
    }

    virtual void updateBoundingEnvelope() = 0;
    size_t directionOfMaxVariance() const;
    size_t findSplitIndex();

    virtual std::pair<SsNode*, SsNode*>  insert(const Point& point) = 0;
    virtual std::pair<SsNode*, SsNode*>  insert(const Point& point, std::string path) = 0;

    bool test(bool isRoot = false) const;
    void print(size_t indent = 0) const;

    //virtual void FNDFTrav(const Point& q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator>& L, NType& Dk) const = 0;

    virtual void kNN(const Point& query, const int& k, pq_type & res, NType& radius) = 0;

    virtual void saveToStream(std::ostream &out) const = 0;
    virtual void loadFromStream(std::istream &in, SsNode* parent) = 0;
};

class SsInnerNode : public SsNode {
private:
    std::vector<Point> getEntriesCentroids() const override;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    SsInnerNode(size_t D) : SsNode(D) {}
//    void draw(sf::RenderWindow& window) const;
//    SsInnerNode() = default;
    int count() const override {
        int c = 0;
        for (auto& child : children) {
            if (child != nullptr) {
                c += child->count();
            }
        }
        return c;
    }
    std::pair<SsNode*, SsNode*> split() override;
    std::vector<SsNode*> children;
    std::pair<std::vector<SsNode *>, std::vector<SsNode *>> splitMeans();
    SsNode* findClosestChild(const Point& target) const;
    bool isLeaf() const override { return false; }
    void updateBoundingEnvelope() override;

    void kNN(const Point& query, const int& k, pq_type & res, NType& radius) override;
    std::pair<SsNode*, SsNode*>  insert(const Point& point) override;
    std::pair<SsNode*, SsNode*>  insert(const Point& point, std::string path) override;


    //void FNDFTrav(const Point& q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator>& L, NType& Dk) const override;

    virtual void saveToStream(std::ostream &out) const override;
    virtual void loadFromStream(std::istream &in, SsNode* parent) override;
};

class SsLeaf : public SsNode {
private:
    std::vector<std::string> paths;

    std::vector<Point> getEntriesCentroids() const override;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    SsLeaf(size_t D) : SsNode(D) {}
//    SsLeaf() = default;
//    void draw(sf::RenderWindow& window) const;
    int count() const override {
        return points.size();
    }
    std::pair<SsNode*, SsNode*> split() override;
    std::vector<Point> points;
    std::pair<std::vector<Point>, std::vector<Point>> splitMeans();
    bool isLeaf() const override { return true; }
    void updateBoundingEnvelope() override;

    std::pair<SsNode*, SsNode*>  insert(const Point& point) override;
    std::pair<SsNode*, SsNode*>  insert(const Point& point, std::string path) override;


    void kNN(const Point& query, const int& k, pq_type & res, NType& radius) override;
    // void FNDFTrav(const Point& q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator>& L, NType& Dk) const override;

    virtual void saveToStream(std::ostream &out) const override;
    virtual void loadFromStream(std::istream &in, SsNode* parent) override;
};

class SsTree {
private:
    SsNode* root;
    SsNode* search(SsNode* node, const Point& target);
    SsNode* searchParentLeaf(SsNode* node, const Point& target);

public:
    size_t D;
    SsTree(size_t D) : root(nullptr), D(D) {}
//    SsTree() : root(nullptr) {}
    ~SsTree() {
        delete root;
    }
    void insert(const Point& point);
    void insert(const Point& point, const std::string& path);
    void build (const std::vector<Point>& points);
    std::vector<std::string> kNNQuery(const Point& center, size_t k) const;

    void print() const;
    void test() const;

    void saveToFile(const std::string &filename) const;
    void loadFromFile(const std::string &filename);
    void printRadius();
    void printChildRadius();
    int count() {
        return root->count();
    }
//    void draw(sf::RenderWindow& window) const;
};

#endif // !SSTREE_H