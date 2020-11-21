#include "MathUtils.h"
#include "gnuplot.h"
#include "DataType.h"

typedef Algebric::MultiDimPoint<double> Vec3;
typedef std::function<Vec3(Vec3,Vec3,double)> Field;

const double G = 1;



class Body {
protected:
    Vec3 position;
    Vec3 velocity;
    double mass;

public:
    Body(double _mass, const Vec3 &pos, const Vec3 &vel) {
        mass = _mass;
        position = pos;
        velocity = vel;
    }
    Body(const Body& b){
        mass = b.mass;
        position =b.position;
        velocity =b.velocity;
    }

    Vec3 &acc_Position() {
        return position;
    }

    Vec3 &acc_getVelocity() {
        return velocity;
    }

    double &acc_Mass() {
        return mass;
    }
};

class World {

    class PhysicalBody : private Body{
    private:
        Field exerted = [](const Vec3& vel,const Vec3& pos,double time)->Vec3{return {0,0,0};}; // by default Nothing!
        Field generated = [this](const Vec3& vel,const Vec3& pos,double time)->Vec3{    // only gravity
            Vec3 relative_pos = pos - position;  // R = r - r'
            return (-G*mass/std::pow(relative_pos.Norm(),3))*relative_pos;
        };
    public:
        PhysicalBody(const Body &b1) : Body(b1) {

        };
    };

private:
    std::vector<PhysicalBody> Bodies;
    double tp;
public:
    World(){

    }

    void AddBody(const Body &b){
        Bodies.emplace_back(b);
    }

};


int main() {
    Body b1(1, {1, 0, 0}, {0, 0, 0});
    Body b2(1, {0, 1, 0}, {0, 0, 0});
    World w;
    w.AddBody(b1);

    return 0;
}