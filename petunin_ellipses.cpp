#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <limits>
#include <iomanip>
#include <sstream>

float M_PI = 3.14159265359;

using namespace std;

struct Point {
    double x, y;
    int label;
};

struct Line {
    double A, B, C;
};

struct Scaling_Factors {
    double scale_x;
    double scale_y;
    double center_x;
    double center_y;
};

pair<vector<Point>, vector<array<double, 3>>> read_Data_From_File(const string& filename) {
    vector<Point> points;
    vector<array<double, 3>> inequalities;
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "Помилка: не вдалося вiдкрити файл " << filename << endl;
        return { points, inequalities };
    }

    bool readingPoints = true;
    while (getline(file, line)) {
        if (line.empty()) continue;

        if (line == "INEQUALITIES:") {
            readingPoints = false;
            continue;
        }

        istringstream iss(line);
        if (readingPoints) {
            double x, y;
            int label;
            if (iss >> x >> y >> label) {
                points.push_back({ x, y, label });
            }
            else {
                cerr << "Помилка: некоректний формат даних для точки у файлi " << filename << endl;
            }
        }
        else {
            double a, b, c;
            if (iss >> a >> b >> c) {
                inequalities.push_back({ a, b, c });
            }
            else {
                cerr << "Помилка: некоректний формат даних для нерiвностi у файлi " << filename << endl;
            }
        }
    }

    file.close();
    return { points, inequalities };
}

pair<Point, Point> estimate_Bounds(const vector<array<double, 3>>& inequalities) {
    double min_x = -1000, max_x = 1000, min_y = -1000, max_y = 1000;

    for (const auto& ineq : inequalities) {
        double a = ineq[0], b = ineq[1], c = ineq[2];

        if (abs(a) > 1e-9 && abs(b) < 1e-9) {
            if (a > 0) max_x = min(max_x, c / a);
            else       min_x = max(min_x, c / a);
        }
        else if (abs(b) > 1e-9 && abs(a) < 1e-9) {
            if (b > 0) max_y = min(max_y, c / b);
            else       min_y = max(min_y, c / b);
        }
    }

    return { {min_x, min_y}, {max_x, max_y} };
}

bool is_Inside_Polygon(const Point& p, const vector<array<double, 3>>& ineq) {
    for (const auto& eq : ineq) {
        if (eq[0] * p.x + eq[1] * p.y > eq[2] + 1e-9) return false;
    }
    return true;
}

vector<Point> generate_Points_In_Polygon(int numPoints, unsigned seed, const vector<array<double, 3>>& inequalities) {
    vector<Point> points;
    mt19937 gen(seed);

    auto [minPt, maxPt] = estimate_Bounds(inequalities);
    uniform_real_distribution<double> x_dist(minPt.x, maxPt.x);
    uniform_real_distribution<double> y_dist(minPt.y, maxPt.y);

    while (points.size() < numPoints) {
        Point p{ x_dist(gen), y_dist(gen), 1 };
        if (is_Inside_Polygon(p, inequalities)) {
            points.push_back(p);
        }
    }
    return points;
}

Line find_Line_Equation(const Point& p1, const Point& p2) {
    return { p2.y - p1.y, p1.x - p2.x, p2.x * p1.y - p1.x * p2.y };
}

Point find_Center(const vector<Point>& points) {
    Point center{ 0, 0 };
    for (const auto& p : points) {
        center.x += p.x;
        center.y += p.y;
    }
    center.x /= points.size();
    center.y /= points.size();
    return center;
}

void rotate_Points(vector<Point>& points, double angle, double center_x, double center_y) {
    double rad = angle * M_PI / 180.0;
    double cos_a = cos(rad);
    double sin_a = sin(rad);

    for (auto& p : points) {
        double x = p.x - center_x;
        double y = p.y - center_y;
        p.x = center_x + x * cos_a - y * sin_a;
        p.y = center_y + x * sin_a + y * cos_a;
    }
}

void rotate_Point(Point& p, double angle, double center_x, double center_y) {
    double rad = angle * M_PI / 180.0;
    double cos_a = cos(rad);
    double sin_a = sin(rad);
    double x = p.x - center_x;
    double y = p.y - center_y;
    p.x = center_x + x * cos_a - y * sin_a;
    p.y = center_y + x * sin_a + y * cos_a;
}

pair<Point, Point> find_Furthest_Points(const vector<Point>& points) {
    double maxDist = 0;
    Point p1 = points[0], p2 = points[1];
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            double dx = points[i].x - points[j].x;
            double dy = points[i].y - points[j].y;
            double dist = dx * dx + dy * dy;
            if (dist > maxDist) {
                maxDist = dist;
                p1 = points[i];
                p2 = points[j];
            }
        }
    }
    return { p1, p2 };
}

double distance_To_Line(const Point& p, const Line& line) {
    return abs(line.A * p.x + line.B * p.y + line.C) / sqrt(line.A * line.A + line.B * line.B);
}

Point find_Projection(const Point& p, const Line& line) {
    double denominator = line.A * line.A + line.B * line.B;
    if (abs(denominator) < 1e-9) return p;
    double x = (line.B * (line.B * p.x - line.A * p.y) - line.A * line.C) / denominator;
    double y = (line.A * (-line.B * p.x + line.A * p.y) - line.B * line.C) / denominator;
    return { x, y };
}

vector<pair<Point, Point>> find_Perpendiculars(const vector<Point>& points, const Line& line) {
    vector<Point> upper, lower;
    for (const auto& p : points) {
        double side = line.A * p.x + line.B * p.y + line.C;
        (side > 0) ? upper.push_back(p) : lower.push_back(p);
    }

    auto findFurthest = [&line](const vector<Point>& group) {
        return *max_element(group.begin(), group.end(),
            [&line](const Point& a, const Point& b) {
                return distance_To_Line(a, line) < distance_To_Line(b, line);
            });
        };

    vector<pair<Point, Point>> result;
    if (!upper.empty()) result.emplace_back(findFurthest(upper), find_Projection(findFurthest(upper), line));
    if (!lower.empty()) result.emplace_back(findFurthest(lower), find_Projection(findFurthest(lower), line));
    return result;
}

Scaling_Factors scale_To_Square(vector<Point>& set1, vector<Point>& set2,
    pair<Point, Point>& mainLine,
    vector<pair<Point, Point>>& perpendiculars) {
    vector<Point> all_points = set1;
    all_points.insert(all_points.end(), set2.begin(), set2.end());
    all_points.push_back(mainLine.first);
    all_points.push_back(mainLine.second);
    for (const auto& p : perpendiculars) {
        all_points.push_back(p.first);
        all_points.push_back(p.second);
    }

    auto [min_x, max_x] = minmax_element(all_points.begin(), all_points.end(),
        [](const Point& a, const Point& b) { return a.x < b.x; });
    auto [min_y, max_y] = minmax_element(all_points.begin(), all_points.end(),
        [](const Point& a, const Point& b) { return a.y < b.y; });

    double width = max_x->x - min_x->x;
    double height = max_y->y - min_y->y;
    double center_x = (min_x->x + max_x->x) / 2;
    double center_y = (min_y->y + max_y->y) / 2;

    double target_size = max(width, height);
    if (target_size < 1e-9) target_size = 1.0;

    double scale_x = (width == 0) ? 1.0 : target_size / width;
    double scale_y = (height == 0) ? 1.0 : target_size / height;

    auto scalePoint = [&](Point& p) {
        p.x = center_x + (p.x - center_x) * scale_x;
        p.y = center_y + (p.y - center_y) * scale_y;
        };

    for (auto& p : set1) scalePoint(p);
    for (auto& p : set2) scalePoint(p);
    scalePoint(mainLine.first);
    scalePoint(mainLine.second);
    for (auto& [p, proj] : perpendiculars) {
        scalePoint(p);
        scalePoint(proj);
    }

    return { scale_x, scale_y, center_x, center_y };
}

vector<Point> compute_Polygon_Vertices(const vector<array<double, 3>>& inequalities) {
    vector<Point> raw_vertices;
    for (size_t i = 0; i < inequalities.size(); ++i) {
        for (size_t j = i + 1; j < inequalities.size(); ++j) {
            auto& a = inequalities[i];
            auto& b = inequalities[j];
            double det = a[0] * b[1] - b[0] * a[1];
            if (abs(det) < 1e-9) continue;

            double x = (b[1] * a[2] - a[1] * b[2]) / det;
            double y = (a[0] * b[2] - b[0] * a[2]) / det;
            Point p{ x, y };
            if (is_Inside_Polygon(p, inequalities)) {
                raw_vertices.push_back(p);
            }
        }
    }
    if (!raw_vertices.empty()) {
        Point center = find_Center(raw_vertices);
        sort(raw_vertices.begin(), raw_vertices.end(), [&](const Point& a, const Point& b) {
            return atan2(a.y - center.y, a.x - center.x) < atan2(b.y - center.y, b.x - center.x);
            });
    }
    return raw_vertices;
}

void write_Data_To_File(const vector<Point>& original_set1,
    const vector<Point>& scaled_set1,
    const vector<Point>& original_set2,
    const vector<Point>& scaled_set2,
    const pair<Point, Point>& originalMainLine,
    const pair<Point, Point>& scaledMainLine,
    const vector<pair<Point, Point>>& originalPerpendiculars,
    const vector<pair<Point, Point>>& scaledPerpendiculars,
    const Point& scaledCenter,
    const Point& originalCenter,
    const Scaling_Factors& factors,
    const vector<Point>& scaledPolygon,
    const vector<Point>& originalPolygon) {
    ofstream out("output_data.txt");
    out << fixed << setprecision(15);

    out << "ScalingFactors:\n" << factors.scale_x << " " << factors.scale_y << "\n\n";
    out << "SquareCenter:\n" << scaledCenter.x << " " << scaledCenter.y << "\n\n";
    out << "OriginalCenter:\n" << originalCenter.x << " " << originalCenter.y << "\n\n";

    out << "CircleRadii:\n";
    for (const auto& p : scaled_set1) {
        double radius = sqrt(pow(p.x - scaledCenter.x, 2) + pow(p.y - scaledCenter.y, 2));
        out << radius << "\n";
    }
    out << "\nEllipseRadii:\n";
    for (const auto& p : scaled_set1) {
        double radius = sqrt(pow(p.x - scaledCenter.x, 2) + pow(p.y - scaledCenter.y, 2));
        out << radius << "\n";
    }

    out << "\nOriginalSet1:\n";
    for (const auto& p : original_set1) out << p.x << " " << p.y << "\n";

    out << "\nScaledSet1:\n";
    for (const auto& p : scaled_set1) out << p.x << " " << p.y << "\n";

    out << "\nOriginalSet2:\n";
    for (const auto& p : original_set2) out << p.x << " " << p.y << "\n";

    out << "\nSet2:\n";
    for (const auto& p : scaled_set2) out << p.x << " " << p.y << "\n";

    out << "\nOriginalMainLine:\n";
    out << originalMainLine.first.x << " " << originalMainLine.first.y << "\n";
    out << originalMainLine.second.x << " " << originalMainLine.second.y << "\n";

    out << "\nMainLine:\n";
    out << scaledMainLine.first.x << " " << scaledMainLine.first.y << "\n";
    out << scaledMainLine.second.x << " " << scaledMainLine.second.y << "\n";

    out << "\nOriginalPerpendiculars:\n";
    for (const auto& [p, proj] : originalPerpendiculars) {
        out << p.x << " " << p.y << " " << proj.x << " " << proj.y << "\n";
    }

    out << "\nPerpendiculars:\n";
    for (const auto& [p, proj] : scaledPerpendiculars) {
        out << p.x << " " << p.y << " " << proj.x << " " << proj.y << "\n";
    }

    out << "\nScaledPolygon:\n";
    for (const auto& p : scaledPolygon) out << p.x << " " << p.y << "\n";

    out << "\nOriginalPolygon:\n";
    for (const auto& p : originalPolygon) out << p.x << " " << p.y << "\n";
}

void input_Inequalities(vector<array<double, 3>>& inequalities, int num_inequalities) {
    cout << "Введiть нерiвностi у форматi: коефiцiєнт_x коефiцiєнт_y вiльний_член (знак >= за замовчуванням)\n";
    inequalities.resize(num_inequalities);
    for (int i = 0; i < num_inequalities; ++i) {
        cout << "Нерiвнiсть " << i + 1 << ": ";
        cin >> inequalities[i][0] >> inequalities[i][1] >> inequalities[i][2];
    }
}

int main() {
    setlocale(LC_ALL, "UKR");

    int Working = 1;
    int Choice = 0;
    unsigned seed;
    int set1_count, total_points, num_inequalities;
    vector<Point> pointsFromFile, original_set1, original_set2;
    vector<array<double, 3>> inequalities;
    mt19937 gen;
    uniform_int_distribution<int> dist(20, 80);

    cout << "------------------------------------------------------------------------------------------------------------------------" << endl;

    while (Working != 0) {
        cout << "Введiть як ви хочете задати данi: " << endl;
        cout << "1 - Зчитування з файлу" << endl;
        cout << "2 - Випадково згенерувати" << endl;
        cout << "0 - Вихiд" << endl;
        cin >> Choice;

        switch (Choice) {
        case 1: {
            auto [points, ineq] = read_Data_From_File("input_data.txt");
            if (!points.empty() || !ineq.empty()) {
                pointsFromFile = points;
                inequalities = ineq;

                original_set1.clear();
                original_set2.clear();
                for (const auto& p : pointsFromFile) {
                    if (p.label == 1) original_set1.push_back(p);
                    else original_set2.push_back(p);
                }
            }
            break;
        }

        case 2: {
            cout << "Введiть загальну кiлькiсть точок: ";
            cin >> total_points;
            cout << "Введiть кiлькiсть нерiвностей (4 або 5): ";
            cin >> num_inequalities;
            while (num_inequalities != 4 && num_inequalities != 5) {
                cout << "Невiрний вибiр. Введiть 4 або 5: ";
                cin >> num_inequalities;
            }

            input_Inequalities(inequalities, num_inequalities);

            seed = chrono::system_clock::now().time_since_epoch().count();
            gen.seed(seed);

            set1_count = dist(gen);
            if (set1_count > total_points) set1_count = total_points;
            original_set1 = generate_Points_In_Polygon(set1_count, gen(), inequalities);
            original_set2 = generate_Points_In_Polygon(total_points - set1_count, gen(), inequalities);

            for (auto& p : original_set1) p.label = 1;
            for (auto& p : original_set2) p.label = -1;

            pointsFromFile = original_set1;
            pointsFromFile.insert(pointsFromFile.end(), original_set2.begin(), original_set2.end());
            break;
        }

        case 0:
            Working = 0;
            break;

        default:
            cout << "Введiть число, якi написанi в списку" << endl;
        }

        if (Choice == 1 || Choice == 2) {
            auto originalMainLine = find_Furthest_Points(original_set1);

            Line originalLineEq = find_Line_Equation(originalMainLine.first, originalMainLine.second);

            auto originalPerpendiculars = find_Perpendiculars(original_set1, originalLineEq);

            vector<Point> scaled_set1 = original_set1;
            vector<Point> scaled_set2 = original_set2;

            Point center = find_Center(scaled_set1);
            double angle = atan2(originalMainLine.second.y - originalMainLine.first.y,
                originalMainLine.second.x - originalMainLine.first.x);
            double angle_deg = -angle * 180.0 / M_PI;

            rotate_Points(scaled_set1, angle_deg, center.x, center.y);
            rotate_Points(scaled_set2, angle_deg, center.x, center.y);
            auto scaledMainLine = originalMainLine;
            rotate_Point(scaledMainLine.first, angle_deg, center.x, center.y);
            rotate_Point(scaledMainLine.second, angle_deg, center.x, center.y);

            Line mainLine = find_Line_Equation(scaledMainLine.first, scaledMainLine.second);
            auto scaledPerpendiculars = find_Perpendiculars(scaled_set1, mainLine);

            Scaling_Factors factors = scale_To_Square(scaled_set1, scaled_set2, scaledMainLine, scaledPerpendiculars);

            Point scaledCenter = find_Center(scaled_set1);
            Point originalCenter = {
                factors.center_x + (scaledCenter.x - factors.center_x) / factors.scale_x,
                factors.center_y + (scaledCenter.y - factors.center_y) / factors.scale_y
            };

            vector<Point> originalPolygon = compute_Polygon_Vertices(inequalities);
            vector<Point> scaledPolygon = originalPolygon;
            for (auto& p : scaledPolygon) rotate_Point(p, angle_deg, center.x, center.y);
            for (auto& p : scaledPolygon) {
                p.x = factors.center_x + (p.x - factors.center_x) * factors.scale_x;
                p.y = factors.center_y + (p.y - factors.center_y) * factors.scale_y;
            }

            write_Data_To_File(original_set1, scaled_set1, original_set2, scaled_set2,
                originalMainLine, scaledMainLine,
                originalPerpendiculars, scaledPerpendiculars,
                scaledCenter, originalCenter,
                factors, scaledPolygon, originalPolygon);
            cout << "Данi успiшно збережено у output_data.txt\n";
            cout << "------------------------------------------------------------------------------------------------------------------------" << endl;
        }
    }

    system("pause");
    return 0;
}
