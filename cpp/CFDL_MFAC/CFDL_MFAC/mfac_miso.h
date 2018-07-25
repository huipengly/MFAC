/**
 * @brief mfac-MISO方法
 * 
 * @file mfac_miso.h
 * @author huipengly
 * @date 2018-07-25
 * @note 无模型自适应控制，多入单出控制器。此类实现了2入1出。
 */
#pragma once
#include <cmath>
#include <vector>
using std::vector;

class MfacMiso
{
public:
    MfacMiso();
    ~MfacMiso();
    vector<double> out(double yd, double y);		//!< SISO/MISO的控制器输出
protected:
    double eta_;           					//!< η
    double mu_;            					//!< μ
    double rho_;							//!< ρ
    double lambda_;							//!< λ
    double epsilon_;						//!< ε
    double M_;                              //!< todo:搞清楚这是啥？
    double phi1_init_value_;                //!< φ(1)的值
    double phi2_init_value_;                //!< φ(1)的值

    vector<double> y_;				        //!< 实际输出值
    vector<double> u1_;				        //!< 控制器输入
    vector<double> u2_;				        //!< 控制器输入
    vector<double> phi1_;			        //!< PPD φ的估算值
    vector<double> phi2_;			        //!< PPD φ的估算值
};

