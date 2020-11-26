#include "MathUtils.h"
#include "gnuplot.h"
#include "DataType.h"
#include "transform.h"
#include "climits"

typedef Algebric::MultiDimPoint<double> Vec3;
typedef std::function<Vec3(Vec3,Vec3,double)> Field;

const double G = 1;


class Body {
    static int i;
protected:
    int tag;
    Vec3 position;
    Vec3 velocity;
    double mass;

public:
    Body(){};
    Body(double _mass, const Vec3 &pos, const Vec3 &vel) {
        i++;
        tag =i;
        mass = _mass;
        position = pos;
        velocity = vel;
    }
    Body(const Body& b){
        mass = b.mass;
        position =b.position;
        velocity =b.velocity;
        tag = b.tag;
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
    int getTag() const{
        return tag;
    }
};
int Body::i = 0;


class World {
    class PhysicalBody : public Body{
    public:
        PhysicalBody(const Body &b1):Body(b1){ // when its lvalue just copy it

        };
        Field generated = [this](const Vec3& vel,const Vec3& pos,double time,int i=0)->Vec3{   // only gravity
            Vec3 relative_pos = pos - this->position;
            if(relative_pos.Norm()<=1e-18)
                return {0,0,0};
            else{
                 // R = r - r'
                 Vec3 p = position;
                return (-G*mass/std::pow(relative_pos.Norm(),3))*relative_pos;
            }
        };
    };
private:
    std::vector<PhysicalBody> Bodies;
    Field net_field = [](const Vec3& vel,const Vec3 &pos,const double& t)->Vec3{return {0,0,0};};
    double tp;
    std::function<Vec3(double)>(*ODESolver)(const std::function<Vec3(Vec3,Vec3,double)> &F_yp_y_Vec3 ,double t0,Vec3 yp0,Vec3 y0, double sTepsize);
public:
    World(){
        tp =0;
    }

    void AddBody(const Body &b){
        int i = Bodies.size();
        Bodies.emplace_back(b);
        Field Copy = net_field;
        net_field = [Copy,this,i](const Vec3 &vel,const Vec3 &pos,const double t)->Vec3{
            PhysicalBody &b = Bodies.at(i);
            return Copy(vel,pos,t) + b.generated(vel,pos,t);
        };
    }
    Field getField(){
        Bodies.at(0).acc_Mass() =0;
        return net_field;
    }
    void init(std::function<Vec3(double)>(*odeSolver)(const std::function<Vec3(Vec3, Vec3, double)> &F_yp_y_Vec3 , double t0, Vec3 yp0, Vec3 y0, double sTepsize)){
        this->ODESolver = odeSolver;
    };

    void ElapseInTime(double stepsize){
        for(auto& b: Bodies){
            b.acc_Position() = ODESolver(net_field,tp,b.acc_getVelocity(),b.acc_Position(),stepsize)(tp+stepsize);
        }
        tp += stepsize;
    }
    void LongRun(double HowMuch_question_Mark,double stepsize){
        int n = std::abs(HowMuch_question_Mark/stepsize) +1;
        stepsize = HowMuch_question_Mark/n;
        for (int i = 0; i < n; ++i) {
            ElapseInTime(stepsize);
        }
    }
    Body acc_Body(const Body& b){
        for(const auto& b1:Bodies){
            if(b1.getTag() == b.getTag()){
                return b1;
            }
        }
        return Body();
    };
    std::vector<Vec3> takeShot(){
        std::vector<Vec3> ans;
        for(auto &b : Bodies)
            ans.push_back(b.acc_Position());
        return ans;
    }
    std::vector<std::vector<double>> takeShotVector(){
        std::vector<std::vector<double>> ans;
        for(auto &b : Bodies)
            ans.push_back(b.acc_Position().getVec());
        return ans;
    }

};


int main() {

    std::function<Algebric::MultiDimPoint<double>(double)> PhaseCurve = [](double t){
        Body b1(1, {1, 0, 0}, {0.5, 0, 0});
        Body b2(1, {0, 1, 0}, {0, 0, 0});
        World w;
        w.AddBody(b1);
        w.AddBody(b2);
        w.init(Calculus::SingleVar::ODE<Vec3,double>::EulerMethodSecondOrder);
        w.LongRun(t,0.01);
        return Algebric::MultiDimPoint<double>(Utils::ConcateVector(w.takeShotVector()));
    };

    plt::ThreeDim::CurvePlot3D<double,double>(PhaseCurve,0,2,100);





    return 0;
}