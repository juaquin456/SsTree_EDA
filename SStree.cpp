#include "SStree.h"






#include "SStree.h"

#define NORMAL

Point calculateMidPoint(std::vector<Point>& points ) {
    Point centroid = points.at(0);
    for (int i=0; i < points.at(0).dim(); i++) {
        NType max = NType::min_value();
        NType min = NType::max_value();
        for (auto p: points) {
            if (p[i] > max) {
                max = p[i];
            }
            if (p[i] < min) {
                min = p[i];
            }
        }
        centroid[i] = (max + min) / 2.0;
    }
    return centroid;
}

Point CalculateMean(std::vector<SsNode*> points) {
    const size_t d = points.at(0)->centroid.dim();
    Point centroid = points.at(0)->centroid.dim();
    for (int i=0; i < d; i++) {
        NType mean = 0;
        for (auto p: points) {
            mean += p->centroid[i];
        }
        mean /= points.size();
        centroid[i] = mean;
    }
    return centroid;
}

Point CalculateMean(std::vector<Point> points) {
    const size_t d = points.at(0).dim();
    Point centroid = points.at(0);
    for (int i=0; i < d; i++) {
        NType mean = 0;
        for (auto p: points) {
            mean += p[i];
        }
        mean /= points.size();
        centroid[i] = mean;
    }
    return centroid;
}



size_t SsNode::directionOfMaxVariance() const {
    NType maxVariance = 0;
    auto directionIndex = 0;
    auto c = this->getEntriesCentroids();
    for (int i = 0; i < this->centroid.dim(); i++) {
        auto t = varianceAlongDirection(c, i);
        if (t > maxVariance) {
            maxVariance = t;
            directionIndex = i;
        }
    }
    return directionIndex;
}

size_t SsNode::findSplitIndex() {
    auto t = this->directionOfMaxVariance();
    this->sortEntriesByCoordinate(t);
    return minVarianceSplit(t);
}

size_t SsNode::minVarianceSplit(size_t coordinateIndex) {
    auto values = this->getEntriesCentroids();
    auto minVariance = NType::max_value();
    auto split_index = Settings::m;
    for (int i = Settings::m; i <= this->getEntriesCentroids().size() - Settings::m; i++) {
        auto var1 = varianceAlongDirection(std::vector<Point>(values.begin(), values.begin() + i), coordinateIndex);
        auto var2 = varianceAlongDirection(std::vector<Point>(values.begin() + i, values.end()), coordinateIndex);
        if (var1 + var2 < minVariance) {
            minVariance = var1 + var2;
            split_index = i;
        }
    }
    return split_index;
}

NType SsNode::varianceAlongDirection(const std::vector<Point> &centroids, size_t direction) const {
    NType mean(0);
    for (auto& c: centroids) {
        mean += c[direction];
    }
    mean /= centroids.size();

    NType variance(0);
    for (auto& c: centroids) {
        variance += pow(c[direction] - mean, 2);
    }
    variance /= centroids.size();
    return variance;
}


void SsTree::build(const std::vector<Point> &points) {


}

SsNode *SsTree::searchParentLeaf(SsNode *node, const Point &target) {
    if (node->isLeaf()) {
        return node;
    }
    else {
        return this->searchParentLeaf(((SsInnerNode*)node)->findClosestChild(target), target);
    }
}

void SsTree::insert(const Point &point) {
    if (this->root == nullptr) {
        this->root = new SsLeaf(D);
        this->root->centroid = point;
        this->root->radius = 0;
    }
    auto p = this->root->insert(point);
    if (p.first != nullptr) {
        this->root = new SsInnerNode(D);
        p.first->parent = this->root;
        p.second->parent = this->root;
        ((SsInnerNode*)root)->children.push_back(p.first);
        ((SsInnerNode*)root)->children.push_back(p.second);
        this->root->updateBoundingEnvelope();
    }
    std::cout << this->root->radius << std::endl;
}

std::vector<std::string> SsTree::kNNQuery(const Point &center, size_t k) const {
    std::priority_queue<std::pair<NType,std::pair<Point, std::string>>, std::vector<std::pair<NType,std::pair<Point, std::string>>>, decltype(comp)> res(comp);
    for (int i = 0; i < k; i++) {
        res.emplace(NType::max_value(), std::make_pair(Point(this->root->centroid), ""));
    }
    auto radius = NType ::max_value();
    this->root->kNN(center, k, res, radius);

    std::vector<std::string> r;

    for (int i=0; i< k; i++) {
        std::cout << res.top().first << "\n";
        r.push_back(res.top().second.second);
        res.pop();
    }

    return r;
}

SsNode *SsInnerNode::findClosestChild(const Point &target) const {
    auto min_dist = NType::max_value();
    auto count_child = NType ::max_value();
    SsNode* res = nullptr;
    for (auto& child: this->children) {
        const auto d = max(distance(child->centroid, target) - child->radius, NType (0));
        if ( d <= min_dist) {
            if (d == 0 and min_dist == 0) {
                if (child->centroid[0] < count_child) {
                    count_child = child->centroid[0];
                    res = child;
                }
            }
            else {
                min_dist = d;
                res = child;
            }
        }
    }

    return res;
}

std::pair<SsNode *, SsNode *> SsInnerNode::insert(const Point &point) {
    auto t = this->findClosestChild(point);
    auto p = t->insert(point);
    if (p.first == nullptr) {
        this->updateBoundingEnvelope();
        return {};
    }
    else {
        this->children.erase(std::remove(this->children.begin(), this->children.end(), t));
        this->children.push_back(p.first);
        this->children.push_back(p.second);

        this->updateBoundingEnvelope();
        if (this->children.size() <= Settings::M){
            return {};
        }
    }
    return this->split();
}

std::pair<SsNode *, SsNode *> SsInnerNode::split() {
#ifdef NORMAL
    auto spl = this->findSplitIndex();
    auto node1 = new SsInnerNode(D);
    node1->children = std::move(std::vector<SsNode*>(this->children.begin(), this->children.begin() + spl));
    node1->updateBoundingEnvelope();
    node1->parent = this->parent;

    auto node2 = new SsInnerNode(D);
    node2->children = std::move(std::vector<SsNode*>(this->children.begin()+spl, this->children.end()));
    node2->updateBoundingEnvelope();
    node2->parent = this->parent;
    return std::make_pair(node1, node2);
#else
    auto r = this->splitMeans();
    if (r.first.empty()) {
        auto spl = this->findSplitIndex();
        auto node1 = new SsInnerNode();
        node1->children = std::vector<SsNode*>(this->children.begin(), this->children.begin() + spl);
        node1->updateBoundingEnvelope();
        node1->parent = this->parent;

        auto node2 = new SsInnerNode();
        node2->children = std::vector<SsNode*>(this->children.begin()+spl, this->children.end());
        node2->updateBoundingEnvelope();
        node2->parent = this->parent;
        return std::make_pair(node1, node2);
    }
    else {
        auto node1 = new SsInnerNode();
        node1->children = r.first;
        node1->updateBoundingEnvelope();
        node1->parent = this->parent;

        auto node2 = new SsInnerNode();
        node2->children = r.second;
        node2->updateBoundingEnvelope();
        node2->parent = this->parent;

        return std::make_pair(node1, node2);
    }

#endif
}

Point calculateMidPoint(std::vector<SsNode*> spheres) {
    Point centroid = spheres.at(0)->centroid;
    for (int i=0; i < spheres.at(0)->centroid.dim(); i++) {
        NType max = NType::min_value();
        NType min = NType::max_value();
        for (auto p: spheres) {
            if (p->centroid[i]+p->radius > max) {
                max = p->centroid[i];
            }
            if (p->centroid[i]-p->radius < min) {
                min = p->centroid[i];
            }
        }
        centroid[i] = (max + min) / 2.0;
    }
    return centroid;
}

void SsInnerNode::updateBoundingEnvelope() {
    auto t = this->getEntriesCentroids();
    if (t.empty()) {
        return;
    }
#ifdef NORMAL
    const auto p = t.at(0);
    std::vector<Point> T;
    NType max_dist = NType ::min_value();
    std::vector<int> support(2);
    for(int i=0; i<t.size(); i++) {
        for(int j=i; j<t.size(); j++) {
            const auto d = distance(t.at(i), t.at(j)) + this->children.at(i)->radius + this->children.at(j)->radius;
            if (d > max_dist or t.size() == 1) {
                max_dist = d;
                support.clear();
                support.push_back(i);
                support.push_back(j);
            }
        }
    }
//    std::cout << support.at(0) << "\t" << support.at(1) << "\t"<< t.size()<< std::endl;
    if (support.at(0) == support.at(1)) {
        this->centroid = t.at(support.at(0));
        this->radius = this->children.at(support.at(0))->radius;
    }
    else {
        const auto diff = t.at(support.at(0)) - t.at(support.at(1));
        this->centroid = (t.at(support.at(0)) + t.at(support.at(1)) + (diff/diff.norm())*(this->children.at(support.at(0))->radius - this->children.at(support.at(1))->radius)) / 2.0; // (A+B + (ra - rb)*|A-B|) / 2
        this->radius = distance(this->centroid, t.at(support.at(0))) + this->children.at(support.at(0))->radius;
    }

    std::priority_queue<std::pair<NType, int>> uncovered;
    for (int i = 0; i < t.size(); i++) {
        auto tmp = distance(this->centroid, t.at(i)) + this->children.at(i)->radius;
        if (!(tmp <= this->radius)) {
            uncovered.emplace(tmp, i);
        }
    }
    while (!uncovered.empty()) {
        auto& tmp = uncovered.top();
        uncovered.pop();
        const auto d = distance(this->centroid, t.at(tmp.second)) + this->children.at(tmp.second)->radius;
        if (!(d <= this->radius)) {
            this->radius = this->radius + 0.01;
            this->centroid = this->centroid + (t.at(tmp.second) - this->centroid).unit()/8;

            while (!uncovered.empty()) {
                uncovered.pop();
            }

            for (int i = 0; i < t.size(); i++) {
                // clean the queue
                auto tmp1 = distance(this->centroid, t.at(i)) + this->children.at(i)->radius;
                if ( tmp1 > this->radius) {
                    uncovered.emplace(tmp1, i);
                }
            }
        }
    }

    auto mn = NType ::max_value();
    for(int i = 0; i<t.size(); i++) {
        const auto d = distance(this->centroid, t.at(i)) + this->children.at(i)->radius;
//        if (d > this->radius) {
//            this->radius = d;
//        }
        if (d < mn)  {
            mn = d;
        }
    }
    this->r_min = mn;
#else
    this->centroid = t[0];
    for (int i=0; i < this->centroid.dim(); i++) {
        NType mean = 0;
        for (auto p: t) {
            mean += p[i];
        }
        mean /= t.size();
        this->centroid[i] = mean;
    }
    auto mx = NType ::min_value();
    auto mn = NType ::max_value();
    for (int i = 0; i < t.size(); i++) {
        auto d = distance(this->centroid, t[i]) + this->children[i]->radius;
        if (d > mx) {
            mx = d;
        }
        if (d < mn) {
            mn = d;
        }
    }
    this->radius = mx;
    this->r_min = mn;

#endif
}

std::vector<Point> SsInnerNode::getEntriesCentroids() const {
    std::vector<Point> p;
    for (auto c: this->children) {
        p.push_back(c->centroid);
    }
    return p;
}

void SsInnerNode::sortEntriesByCoordinate(size_t coordinateIndex) {
    std::sort(this->children.begin(), this->children.end(), [coordinateIndex](SsNode* &e1, SsNode* &e2){return e1->centroid[coordinateIndex] < e2->centroid[coordinateIndex];});
}

void
SsInnerNode::kNN(const Point &query, const int &k, pq_type &res, NType &radius) {
    for (auto c: children) {
        auto d_tf = distance(query, c->centroid);
        auto d = d_tf - c->radius;
        if (d < 0) d = 0;
        if (d < radius and !(radius + c->radius < d_tf) and !(radius + d_tf < c->r_min)) {
            c->kNN(query, k, res, radius);
        }
    }
}

int no_reinsert = 0;

std::pair<SsNode *, SsNode *> SsLeaf::insert(const Point &point) {
    if (std::find(this->points.begin(), this->points.end(),point) != this->points.end()) {
        return {};
    }
//    if (this->points.size() == Settings::M and no_reinsert and this->parent != nullptr  ) {
//        // reinsert
//        std::vector<Point> to_reinsert(this->points.begin(), this->points.end());
//        Point mid = calculateMidPoint(this->points);
//        std::sort(to_reinsert.begin(), to_reinsert.end(), [mid](Point &e1, Point &e2){return distance(e1, mid) < distance(e2, mid);});
//        this->points.clear();
//        this->points = std::vector<Point>(to_reinsert.begin(), to_reinsert.begin() + Settings::m);
//        this->updateBoundingEnvelope();
//        // iter the rest
//        no_reinsert = false;
//        for (int i = Settings::m; i < to_reinsert.size(); i++) {
//            this->parent->insert(to_reinsert[i]);
//        }
//        no_reinsert = true;
//        std::cout << "reinserted\n";
//    }
    this->points.push_back(point);
    this->updateBoundingEnvelope();
    if (this->points.size() <= Settings::M) {
        return {};
    }
    return this->split();
}

std::pair<SsNode *, SsNode *> SsLeaf::split() {
#ifdef NORMAL_LEAF

    auto spl = this->findSplitIndex();
    auto node1 = new SsLeaf();
    node1->points = std::vector<Point>(this->points.begin(), this->points.begin() + spl);
    node1->updateBoundingEnvelope();
    node1->parent = this->parent;

    auto node2 = new SsLeaf();
    node2->points = std::vector<Point>(this->points.begin()+spl, this->points.end());
    node2->updateBoundingEnvelope();
    node2->parent = this->parent;
    return std::make_pair(node1, node2);


#else
//    auto r = this->splitMeans();
//    if (r.first.empty()) {
    auto spl = this->findSplitIndex();
    auto node1 = new SsLeaf(D);
    node1->points = std::move(std::vector<Point>(this->points.begin(), this->points.begin() + spl));
    if (!paths.empty())
        node1->paths = std::move(std::vector<std::string>(this->paths.begin(), this->paths.begin() + spl));
    node1->updateBoundingEnvelope();
    node1->parent = this->parent;

    auto node2 = new SsLeaf(D);
    node2->points = std::move(std::vector<Point>(this->points.begin()+spl, this->points.end()));
    if (!paths.empty())
        node2->paths = std::move(std::vector<std::string>(this->paths.begin() + spl, this->paths.end()));
    node2->updateBoundingEnvelope();
    node2->parent = this->parent;
    return std::make_pair(node1, node2);
//    }
//    else {
//        auto node1 = new SsLeaf();
//        node1->points = r.first;
//        node1->updateBoundingEnvelope();
//        node1->parent = this->parent;
//
//        auto node2 = new SsLeaf();
//        node2->points = r.second;
//        node2->updateBoundingEnvelope();
//        node2->parent = this->parent;
//
//        return std::make_pair(node1, node2);
//    }
#endif
}

std::vector<Point> SsLeaf::getEntriesCentroids() const {
    return this->points;
}


void SsLeaf::updateBoundingEnvelope() {
    auto t = this->getEntriesCentroids();
    if (t.empty()) {
        return;
    }

#ifdef NORMAL
    const auto p = t.at(0);
    std::vector<Point> T;
    NType max_dist = NType ::min_value();
    std::vector<int> support(2);
    for(int i=0; i<t.size(); i++) {
        for(int j=i; j<t.size(); j++) {
            const auto d = distance(t.at(i), t.at(j));
            if (d > max_dist or t.size() == 1) {
                max_dist = d;
                support.clear();
                support.push_back(i);
                support.push_back(j);
            }
        }
    }
    this->centroid = (t.at(support.at(0)) + t.at(support.at(1))) / 2.0;
    this->radius = distance(this->centroid, t.at(support.at(1)));
    std::priority_queue<std::pair<NType, int>> uncovered;
    for (int i = 0; i < t.size(); i++) {
        auto tmp = distance(this->centroid, t.at(i));
        if (!(tmp <= this->radius)) {
            uncovered.emplace(tmp, i);
        }
    }

//    if (t.size() == 3) {
//        std::cout << "here\n";
//        std::cout << uncovered.size() << "\n";
//    }
    while (!uncovered.empty()) {
        auto& tmp = uncovered.top();
        uncovered.pop();
        const auto d = distance(this->centroid, t.at(tmp.second));
        if (!(d <= this->radius)) {
            this->radius = this->radius + 0.01;
            this->centroid = this->centroid + (t.at(tmp.second) - this->centroid).unit()/8;

            while (!uncovered.empty()) {
                uncovered.pop();
            }
            for (int i = 0; i < t.size(); i++) {
                if (auto tmp1 = distance(this->centroid, t.at(i)) > this->radius) {
                    uncovered.emplace(tmp1, i);
                }
            }
        }
    }

    auto mn = NType ::max_value();
    for(const auto & i : t) {
        const auto d = distance(this->centroid, i);
        if (d < mn)  {
            mn = d;
        }
    }
    this->r_min = mn;

#else
    this->centroid = t[0];
    for (int i=0; i < this->centroid.dim(); i++) {
        NType mean = 0;
        for (auto p: t) {
            mean += p[i];
        }
        mean /= t.size();
        this->centroid[i] = mean;
    }
    auto mx = NType ::min_value();
    auto mn = NType ::max_value();
    for (const auto & i : t) {
        auto d = distance(this->centroid, i);
        if (d > mx) {
            mx = d;
        }
        if (d < mn) {
            mn = d;
        }
    }
    this->radius = mx;
    this->r_min = mn;

#endif
}

void SsLeaf::sortEntriesByCoordinate(size_t coordinateIndex) {
    std::vector<std::pair<Point, std::string>> entries;
    for (size_t i = 0; i < this->points.size(); i++) {
        entries.emplace_back(this->points[i], this->paths[i]);
    }
    std::sort(entries.begin(), entries.end(), [coordinateIndex](std::pair<Point, std::string>& a, std::pair<Point, std::string>& b) { return a.first[coordinateIndex] < b.first[coordinateIndex]; });

    this->points = std::vector<Point>(entries.size());
    this->paths = std::vector<std::string>(entries.size());

    for (size_t i = 0; i < entries.size(); i++) {
        this->points[i] = entries[i].first;
        this->paths[i] = entries[i].second;
    }
//std::sort(this->points.begin(), this->points.end(), [coordinateIndex](Point& a, Point& b) { return a[coordinateIndex] < b[coordinateIndex]; });
}


void SsLeaf::kNN(const Point& query, const int &k, pq_type &res, NType& radius) {
    auto d_qmp = distance(this->centroid, query);
    for (int i = 0; i < this->points.size(); i++) {
        auto p = this->points[i];
        auto d = distance(this->centroid, p);
        if (d < radius and !(radius + d < d_qmp) and !(radius + d_qmp < d)) {
            auto d1 = distance(query, p);
            res.emplace(d1, std::make_pair(p, this->paths.at(i)));
            res.pop();
            radius = res.top().first;
        }
    }
}

std::pair<std::vector<Point>, std::vector<Point>> SsLeaf::splitMeans() {
    auto p = this->getEntriesCentroids();
    const int k = 2;
    Point centroids[k] = {p.at(0), p.at(1)};
    // 2 - Means
    std::vector<Point> copy[k];
    std::vector<Point> cluster[k];
    for (int s = 0; s < 500; s++) {
        cluster[0].clear();
        cluster[1].clear();
        for (const auto& point: p) {
            const auto d0 = distance(point, centroids[0]);
            const auto d1 = distance(point, centroids[1]);

            if (d0 < d1) {
                cluster[0].push_back(point);
            }
            else {
                cluster[1].push_back(point);
            }
        }
        if (cluster[0].size() < Settings::m or cluster[1].size() < Settings::m) {
            return std::make_pair(copy[0], copy[1]);
        }
        copy[0] = cluster[0];
        copy[1] = cluster[1];
        // RECALCULATE CENTROID
        centroids[0] = CalculateMean(cluster[0]);
        centroids[1] = CalculateMean(cluster[1]);
    }
    return std::make_pair(cluster[0], cluster[1]);
}

std::pair<SsNode *, SsNode *> SsLeaf::insert(const Point &point, std::string path) {
    if (std::find(this->points.begin(), this->points.end(),point) != this->points.end()) {
        return {};
    }
//    if (this->points.size() == Settings::M and no_reinsert % 50 == 0 and this->parent != nullptr  ) {
//        // reinsert
//        std::vector<Point> to_reinsert(this->points.begin(), this->points.end());
//        std::vector<std::string> path_reinsert(this->paths.begin(), this->paths.end());
//        Point mid = calculateMidPoint(this->points);
//        std::sort(to_reinsert.begin(), to_reinsert.end(), [mid](Point &e1, Point &e2){return distance(e1, mid) < distance(e2, mid);});
//        this->points.clear();
//        this->points = std::vector<Point>(to_reinsert.begin(), to_reinsert.begin() + Settings::m);
//        this->paths = std::vector<std::string>(this->paths.begin(), this->paths.begin() + Settings::m);
//        this->updateBoundingEnvelope();
//        // iter the rest
//        no_reinsert = false;
//        for (int i = Settings::m; i < to_reinsert.size(); i++) {
//            this->parent->insert(to_reinsert[i], path_reinsert[i]);
//        }
//        no_reinsert++;
//        std::cout << "reinserted\n";
//    }
    this->points.push_back(point);
    this->paths.push_back(path);
    this->updateBoundingEnvelope();
    if (this->points.size() <= Settings::M) {
        return {};
    }
    return this->split();
}

const float norm = 1.7;
const float offset_x = 170;
const float offset_y = 250;
//void SsLeaf::draw(sf::RenderWindow &window) const {
//    sf::CircleShape circle(radius.getValue()/norm);
//    circle.setPosition((centroid[0].getValue() - radius.getValue())/norm + offset_x, offset_y + (centroid[1].getValue() - radius.getValue())/norm);
//    circle.setFillColor(sf::Color::Transparent);
//    circle.setOutlineThickness(2);
//    circle.setOutlineColor(sf::Color::Blue);
//    window.draw(circle);
//
//    for (const auto& point : points) {
//        sf::CircleShape circle(2);
//        circle.setPosition(point[0].getValue()/norm + offset_x, point[1].getValue()/norm + offset_y);
//        circle.setFillColor(sf::Color::Green);
//        window.draw(circle);
//    }
//}


std::pair<std::vector<SsNode *>, std::vector<SsNode *>> SsInnerNode::splitMeans() {
    auto p = this->getEntriesCentroids();
    const int k = 2;
    Point centroids[k] = {p.at(0), p.at(1)};
    // 2 - Means

    std::vector<SsNode*> copy[k];
    std::vector<SsNode*> cluster[k];
    // k-means with minimum Settings::m points per cluster
    for (int s = 0; s < 500; s++) {
        cluster[0].clear();
        cluster[1].clear();
        for (int i = 0; i < p.size(); i++) {
            const auto d0 = distance(p.at(i), centroids[0]);
            const auto d1 = distance(p.at(i), centroids[1]);
            if (d0 < d1) {
                cluster[0].push_back(this->children.at(i));
            }
            else {
                cluster[1].push_back(this->children.at(i));
            }
        }
        if (cluster[0].size() < Settings::m or cluster[1].size() < Settings::m) {
            return std::make_pair(copy[0], copy[1]);
        }
        copy[0] = cluster[0];
        copy[1] = cluster[1];
        // RECALCULATE CENTROID
        centroids[0] = CalculateMean(cluster[0]);
        centroids[1] = CalculateMean(cluster[1]);
    }
    return std::make_pair(cluster[0], cluster[1]);
}



std::pair<SsNode *, SsNode *> SsInnerNode::insert(const Point &point, std::string path) {
    auto t = this->findClosestChild(point);
    auto p = t->insert(point, path);
    if (p.first == nullptr) {
        this->updateBoundingEnvelope();
        return {};
    }
    else {
        this->children.erase(std::remove(this->children.begin(), this->children.end(), t));
        this->children.push_back(p.first);
        this->children.push_back(p.second);

        this->updateBoundingEnvelope();
        if (this->children.size() <= Settings::M){
            return {};
        }
    }
    return this->split();
}


//void SsInnerNode::draw(sf::RenderWindow &window) const {
//    sf::CircleShape circle(radius.getValue()/norm);
//    circle.setPosition(offset_x + (centroid[0].getValue() - radius.getValue())/norm, offset_y + (centroid[1].getValue() - radius.getValue())/norm);
//    circle.setFillColor(sf::Color::Transparent);
//    circle.setOutlineThickness(2);
//    circle.setOutlineColor(sf::Color::Red);
//    window.draw(circle);
//
//    for (const auto& child : children) {
//        child->draw(window);
//    }
//}


void SsTree::insert(const Point &point, const std::string &path) {
    if (this->root == nullptr) {
        this->root = new SsLeaf(D);
        this->root->centroid = point;
        this->root->radius = 0;
    }
    auto p = this->root->insert(point, path);
    if (p.first != nullptr) {
        this->root = new SsInnerNode(D);
        p.first->parent = this->root;
        p.second->parent = this->root;
        ((SsInnerNode*)root)->children.push_back(p.first);
        ((SsInnerNode*)root)->children.push_back(p.second);
        this->root->updateBoundingEnvelope();
    }
//    std::cout << this->root->radius << std::endl;
}

void SsTree::printRadius() {
    std::cout << "radio del root :" << this->root->radius << "\n";
}

void SsTree::printChildRadius() {
    for (auto c: ((SsInnerNode*)root)->children) {
        std::cout << "radio del hijo :" << c->radius << "\n";
    }
}
//
//void SsTree::draw(sf::RenderWindow &window) const {
//    root->draw(window);
//}














bool SsNode::test(bool isRoot) const {
    size_t count = 0;
    if (this->isLeaf()) {
        const SsLeaf* leaf = dynamic_cast<const SsLeaf*>(this);
        count = leaf->points.size();

        // Verificar si los puntos están dentro del radio del nodo
        for (const Point& point : leaf->points) {
            if (distance(this->centroid, point) > this->radius) {
                std::cout << "Point outside node radius detected." << std::endl;
                return false;
            }
        }
    } else {
        const SsInnerNode* inner = dynamic_cast<const SsInnerNode*>(this);
        count = inner->children.size();

        // Verificar si los centroides de los hijos están dentro del radio del nodo padre
        for (const SsNode* child : inner->children) {
            if (distance(this->centroid, child->centroid) > this->radius) {
                std::cout << "Child centroid outside parent radius detected." << std::endl;
                return false;
            }
            // Verificar recursivamente cada hijo
            if (!child->test()) {
                return false;
            }
        }
    }

    // Comprobar la validez de la cantidad de hijos/puntos
    if (!isRoot && (count < Settings::m || count > Settings::M)) {
        std::cout << "Invalid number of children/points detected." << std::endl;
        return false;
    }

    // Comprobar punteros de parentezco, salvo para el nodo raíz
    if (!isRoot && !parent) {
        std::cout << "Node without parent detected." << std::endl;
        return false;
    }

    return true;
}

void SsTree::test() const {
    bool result = root->test();

    if (root->parent) {
        std::cout << "Root node parent pointer is not null!" << std::endl;
        result = false;
    }

    if (result) {
        std::cout << "SS-Tree is valid!" << std::endl;
    } else {
        std::cout << "SS-Tree has issues!" << std::endl;
    }
}


void SsNode::print(size_t indent) const {
    for (size_t i = 0; i < indent; ++i) {
        std::cout << "  ";
    }

    // Imprime información del nodo.
    std::cout << "Centroid: " << centroid << ", Radius: " << radius;
    if (isLeaf()) {
        const SsLeaf* leaf = dynamic_cast<const SsLeaf*>(this);
        std::cout << ", Points: [ ";
        for (const Point& p : leaf->points) {
            std::cout << p << " ";
        }
        std::cout << "]";
    } else {
        std::cout << std::endl;
        const SsInnerNode* inner = dynamic_cast<const SsInnerNode*>(this);
        for (const SsNode* child : inner->children) {
            child->print(indent + 1); 
        }
    }
    std::cout << std::endl;
}
void SsTree::print() const {
    if (root) {
        root->print();
    } else {
        std::cout << "Empty tree." << std::endl;
    }
}
void SsLeaf::saveToStream(std::ostream &out) const {
    // Guardar centroid
    centroid.saveToFile(out, D);

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char*>(&radius_), sizeof(radius_));

    // Guardar el numero de puntos
    size_t numPoints = points.size();
    out.write(reinterpret_cast<const char*>(&numPoints), sizeof(numPoints));

    // Guardar los puntos
    for (const auto& point : points) {
        point.saveToFile(out, D);
    }

    // Guardar las rutas (paths)
    size_t numPaths = paths.size();
    out.write(reinterpret_cast<const char*>(&numPaths), sizeof(numPaths));
    for (const auto& p : paths) {
        size_t pathLength = p.size();
        out.write(reinterpret_cast<const char*>(&pathLength), sizeof(pathLength));
        out.write(p.c_str(), (long) pathLength);
    }
}

void SsInnerNode::saveToStream(std::ostream &out) const {
    // Guardar centroid
    centroid.saveToFile(out, D);

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char*>(&radius_), sizeof(radius_));

    // Guardar si apunta a nodos hoja
    bool pointsToLeafs = children[0]->isLeaf();
    out.write(reinterpret_cast<const char*>(&pointsToLeafs), sizeof(pointsToLeafs));

    // Guardar la cantidad de hijos para saber cuántos nodos leer después
    size_t numChildren = children.size();
    out.write(reinterpret_cast<const char*>(&numChildren), sizeof(numChildren));

    // Guardar los hijos
    for (const auto& child : children) {
        child->saveToStream(out);
    }
}

void SsInnerNode::loadFromStream(std::istream &in, SsNode* parent) {
    this->parent = parent;

    // Leer centroid
    centroid.readFromFile(in, D);

    // leer el valor del radio
    float radius_ = 0;
    in.read(reinterpret_cast<char*>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // leer si apunta a hojas o nodos internos
    bool pointsToLeaf = false;
    in.read(reinterpret_cast<char*>(&pointsToLeaf), sizeof(pointsToLeaf));

    // leer cantidad de hijos
    size_t numChildren;
    in.read(reinterpret_cast<char*>(&numChildren), sizeof(numChildren));

    // leer hijos
    for (size_t i = 0; i < numChildren; ++i) {
        SsNode* child = pointsToLeaf ? static_cast<SsNode*>(new SsLeaf(D)) : static_cast<SsNode*>(new SsInnerNode(D));
        child->loadFromStream(in, this);
        children.push_back(child);
    }
}


void SsLeaf::loadFromStream(std::istream &in,  SsNode* parent) {
    this->parent = parent;

    // Leer centroid
    centroid.readFromFile(in, D);

    // Leer radio
    float radius_ = 0;
    in.read(reinterpret_cast<char*>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // Leer numero de puntos
    size_t numPoints;
    in.read(reinterpret_cast<char*>(&numPoints), sizeof(numPoints));

    // Leer puntos
    points.resize(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        points[i].readFromFile(in, D);
    }

    // Leer rutas (paths)
    size_t numPaths;
    in.read(reinterpret_cast<char*>(&numPaths), sizeof(numPaths));
    paths.resize(numPaths);
    for (size_t i = 0; i < numPaths; ++i) {
        size_t pathLength;
        in.read(reinterpret_cast<char*>(&pathLength), sizeof(pathLength));
        char* buffer = new char[pathLength + 1];
        in.read(buffer, (long) pathLength);
        buffer[pathLength] = '\0';
        paths[i] = std::string(buffer);
        delete[] buffer;
    }
}

void SsTree::saveToFile(const std::string &filename) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Cannot open file for writing");
    }

    // Guardar las dimensiones de la estructura
    out.write(reinterpret_cast<const char*>(&D), sizeof(D));

    // Guardar si el root es hija o nodo interno
    bool isLeaf = root->isLeaf();
    out.write(reinterpret_cast<const char*>(&isLeaf), sizeof(isLeaf));

    // Guardar el resto de la estructura
    root->saveToStream(out);
    out.close();
}

void SsTree::loadFromFile(const std::string &filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot open file for reading");
    }
    if (root) {
        delete root;
        root = nullptr;
    }

    // Aquí se asume que el primer valor determina las dimensiones
    in.read(reinterpret_cast<char*>(&D), sizeof(D));

    // El segundo valor determina si el root es hoja
    bool isLeaf;
    in.read(reinterpret_cast<char*>(&isLeaf), sizeof(isLeaf));
    if (isLeaf) {
        root = new SsLeaf(D);
    } else {
        root = new SsInnerNode(D);
    }
    root->loadFromStream(in, nullptr);
    in.close();
}