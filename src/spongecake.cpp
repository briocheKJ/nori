
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/common.h>

NORI_NAMESPACE_BEGIN

class SpongeCake : public BSDF {
public:
    SpongeCake(const PropertyList& propList) {
        m_alpha = propList.getFloat("alpha", 0.1f);

        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        m_kd = propList.getColor("kd", Color3f(0.5f));

        m_ks = 1 - m_kd.maxCoeff();
    }
    float D(const Vector3f& omega) const {
        // 计算 S^{-1} 矩阵
        Eigen::Matrix3f S_inv; //类纤维
        S_inv << 1.0f, 0, 0,
            0, 1.0f, 0,
            0, 0, 1.0f / (m_alpha * m_alpha);

        // 计算 q
        float q = omega.transpose() * S_inv * omega;

        // 计算 S 矩阵
        Eigen::Matrix3f S;
        S << 1.0f, 0, 0,
            0, 1.0f, 0,
            0, 0, m_alpha * m_alpha;

        // 计算 sigma
        float sigma = std::sqrt(omega.transpose() * S * omega);

        // 使用 q 和 sigma 计算 D(omega)
        float D_omega = 1.0f / (M_PI * m_alpha * q * q);

        return D_omega;
    }
    float Lambda(const Vector3f& w) const {
        // TODO: 实现 Lambda 函数
        Vector3f n = Vector3f(0.0f, 0.0f, 1.0f);
        return 1.0f/(w.dot(n) / (w.norm() * n.norm()));  // 临时值
    }
    float G(const Vector3f& wi, const Vector3f& wo) const {
        float T = 1.0f;  // 假设 T 为 1，根据你的情况来调整
        float rho = 1.0f;  // 假设 rho 为 1，根据你的情况来调整
        float expTerm = std::exp(-T * rho * (Lambda(wi) + Lambda(wo)));
        return (1.0f - expTerm) / (Lambda(wi) + Lambda(wo));
    }


    // Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord& bRec) const {

        Vector3f wi = bRec.wi.normalized();
        Vector3f wo = bRec.wo.normalized();

        if (Frame::cosTheta(wo) <= 0.0f)
            return Color3f(0.0f);

        Vector3f wh = (wi + wo).normalized();

        Vector3f n = Vector3f(0.0f, 0.0f, 1.0f);

        auto b = 0.0f, c = 0.0f;
        float cosine_theta_v, tan_theta_v;

        float cosThetai = (wh).dot(wi);
        float fresnelterm = fresnel(cosThetai, m_extIOR, m_intIOR);
        float D = this->D(wh);

        //for incoming direction - masking in
        if ((wi.dot(wh) / wi.dot(n)) > 0.0f)
            c = 1.0f;
        else
            c = 0.0f;

        cosine_theta_v = wi.dot(n) / (wi.norm() * n.norm());
        tan_theta_v = (1 - (cosine_theta_v * cosine_theta_v)) / (cosine_theta_v * cosine_theta_v);
        b = 1 / (m_alpha * tan_theta_v);


        //for outgoing direction - masking out
        if ((wo.dot(wh) / wo.dot(n)) > 0.0f)
            c = 1.f;
        else
            c = 0.f;

        cosine_theta_v = wo.dot(n) / (wo.norm() * n.norm());
        tan_theta_v = (1 - (cosine_theta_v * cosine_theta_v)) / (cosine_theta_v * cosine_theta_v);

        b = 1 / (m_alpha * tan_theta_v);

        float G_term = G(wi, wo);
        auto part1 = m_kd / M_PI;
        auto part2 = (m_ks * D * fresnelterm * G_term) / (4 * Frame::cosTheta(wi) * Frame::cosTheta(wo) * Frame::cosTheta(wh));
        return (part1 + part2);
    }


    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord& bRec) const {
        Vector3f wi = bRec.wi.normalized();
        Vector3f wo = bRec.wo.normalized();

        if (Frame::cosTheta(wo) <= 0.0f)
            return 0.0f;

        Vector3f wh = (wi + wo).normalized();

        float jh = 1.0f / (4.0f * (wh.dot(wo)));
        float D = this->D(wh);
        float pdf = (m_ks * D * jh) + (((1.0f - m_ks) * Frame::cosTheta(wo)) / M_PI);
        return pdf;
    }


    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord& bRec, const Point2f& _sample) const {
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);

        Color3f sampledbsdf;
        if (Frame::cosTheta(bRec.wi) <= 0.0f)
            return Color3f(0.0f);
        auto sample = _sample;
        //check if specular or diffuse
        if (_sample.y() < m_ks)
        {
            //specular

            //scaling and offsetting the sample; range of _sample.y() is [0, m_ks] and want it to [0,1]
            float scaledsample = _sample.y() / m_ks;
            Point2f scaledpoint(_sample.x(), scaledsample);
            //sample a normal using squareToBeckmann as specular
            Vector3f whpoint = Warp::squareToBeckmann(scaledpoint, m_alpha);
            Vector3f wh = whpoint.normalized();

            //now we reflect wi using n in order to generate the outgoing direction wo
            Vector3f wo = -bRec.wi + (2 * (wh.dot(bRec.wi)) * wh);
            bRec.wo = wo;
            if (Frame::cosTheta(wo) <= 0)
                return Color3f(0.f);
            sampledbsdf = (eval(bRec) * Frame::cosTheta(bRec.wo)) / pdf(bRec);
        }
        else
        {
            //diffuse

            //scaling and offsetting the sample; range of _sample.y() is [0, m_ks] and want it to [0,1]
            float scaledsample = (_sample.y() - m_ks) / (1 - m_ks);
            Point2f scaledpoint(_sample.x(), scaledsample);
            //sample a normal using squareToCosineHemisphere as specular
            Vector3f wopoint = Warp::squareToCosineHemisphere(scaledpoint);
            bRec.wo = wopoint;
            if (Frame::cosTheta(wopoint) <= 0)
                return Color3f(0.0f);
            sampledbsdf = (eval(bRec) * Frame::cosTheta(bRec.wo)) / pdf(bRec);
        }
        return sampledbsdf;
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "SpongeCake[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(SpongeCake, "spongecake");
NORI_NAMESPACE_END