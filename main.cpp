#include <iostream>
#include <vector>
#include <random>
#include "SStree.h"

int main() {
    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1000.0, 1000.0);

    std::vector<Point> points(500, Point(50));
    for (int i = 0; i < 500; i++) {
        for (int j = 0; j < 50; j++) {
            points[i][j] = dis(gen);
        }
    }

    // Insert the points into the SStree
    SsTree tree(50);
    for (auto& point : points) {
        tree.insert(point);
    }

    std::string filename = "sstree.dat";
    tree.saveToFile(filename);

    tree = SsTree();    // Clean the tree
    tree.loadFromFile(filename);

    tree.kNNQuery(Point(50), 10);
    return 0;
}

//int main() {
//    SsTree tree(448);
////    tree.D = 448;
//    tree.loadFromFile("embbeding.dat");
//    tree.print();
//    tree.kNNQuery(Point(448), 10);
//}