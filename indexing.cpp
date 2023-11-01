#include <iostream>
#include <vector>
#include <random>
#include "Point.h"
#include "SStree.h"
//#include <hdf5/serial/H5Cpp.h>
#include <nlohmann/json.hpp>
#include <fstream>

// g++ -std=c++11 -fopenmp -o indexing indexing.cpp -lfaiss -lopenblas -L/path/to/faiss/lib -I/path/to/faiss/include -I/path/to/nlohmann

struct ImageData {
    std::vector<Point> embeddings;
    std::vector<std::string> paths;
};

ImageData readEmbeddingsFromJson(const std::string& FILE_NAME) {
    ImageData data;

    try {
        std::ifstream file(FILE_NAME);
        if (!file.is_open()) {
            throw std::runtime_error("Unable to open JSON file.");
        }

        nlohmann::json jsonData;
        file >> jsonData;
        
        std::vector<std::vector<float>> features;
        for (const auto& featureList : jsonData["features"]) {
            std::vector<float> tempFeature = featureList;
            features.push_back(tempFeature);
        }

        data.paths = jsonData.at("paths").get<std::vector<std::string>>();

        if (features.size() != data.paths.size()) {
            throw std::runtime_error("The number of features does not match the number of paths.");
        }

        for (const auto& feature : features) {
            Point embedding(feature.size());
            for (size_t j = 0; j < feature.size(); ++j) {
                embedding[j] = feature[j];
            }
            data.embeddings.push_back(embedding);
        }
        
        file.close();

        for (const auto& path : data.paths) {
            if (path.empty()) {
                throw std::runtime_error("Uno de los paths está vacío.");
            }
//            std::cout << "Guardado path: " << path << std::endl;
        }

    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return data;
}


//int main() {
//    const std::string FILE_NAME("embedding.json");
//    ImageData data = readEmbeddingsFromJson(FILE_NAME);
//
//    SsTree tree(448);
//    for (size_t i = 0; i < data.embeddings.size(); ++i) {
//        tree.insert(data.embeddings[i], data.paths[i]);
//        if (i%1000 == 0) tree.printRadius();
//    }
//    tree.printRadius();
//    tree.printChildRadius();
//    std::string filename = "embbeding.dat";
//    tree.saveToFile(filename);
//    tree.kNNQuery(Point(448), 10);
//    return 0;
//}

int main() {
    const std::string FILE_NAME("embedding.json");
    ImageData data = readEmbeddingsFromJson(FILE_NAME);

    SsTree tree(448);
    tree.loadFromFile("embbeding.dat");
    std::cout << "Info del árbol cargado: " << std::endl;
    tree.printRadius();
    tree.printChildRadius();
    std::cout << "--------------------------" << std::endl;
    std::cout << "10 - nn query with index: point(0, 0, 0....)" << std::endl;
    tree.kNNQuery(Point(448), 10);
    std::cout << "10 - nn query with exhaustive search: point(0, 0, 0....)" << std::endl;
    std::greater<NType> greater;
    std::priority_queue<NType, std::vector<NType > ,decltype(greater)> distances(greater);
    for(size_t i = 0; i < data.embeddings.size(); ++i) {
        distances.push(distance(Point(448), data.embeddings[i]));
    }
    for (size_t i = 0; i < 10; ++i) {
        std::cout << distances.top() << std::endl;
        distances.pop();
    }
    std::cout << "--------------------------" << std::endl;
}
