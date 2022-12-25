#include "distance_calculator.h"

Distance_Calculator::Distance_Calculator()
{
    stat_max = 0;
    stat_avg = 0;
    stat_min = std::numeric_limits<double>::max();
    count_external = 0;
}

Distance_Calculator::~Distance_Calculator()
{

}