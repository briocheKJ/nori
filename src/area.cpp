#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/mesh.h>
#include <Eigen/Geometry>
#include <nori/warp.h>


NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter
{
public:
    AreaLight(const PropertyList& propList)
    {
        m_radiance = propList.getColor("radiance", Color3f(1.0f));
    }

    Color3f eval(const EmitterQueryRecord& eRec) const
    {
        return m_radiance;
    }

    float pdf(Mesh* mesh, const EmitterQueryRecord& eRec) const
    {
        //return 1.0f / mesh->getArea();
        return 1.0 / mesh->getTotalArea();
    }

    /// Draw a sample from the Emitter model


    Point3f sample(Mesh* mesh, EmitterQueryRecord& eRec, Sampler* sampler) const
    {
        //lRec.origin is the point that the ray from eye intersects on the scene; point x on the scene where we want to find the radiance
        EmitterQueryRecord SQR(eRec.lightPos);
        const MatrixXf& m_V = mesh->getVertexPositions();
        const MatrixXf& m_N = mesh->getVertexNormals();
        const MatrixXf& m_UV = mesh->getVertexTexCoords();
        const MatrixXu& m_F = mesh->getIndices();

        const Point2f& sample = sampler->next2D();
        //take random sample of the index of the triangle
        DiscretePDF m_pdf1 = mesh->getTrianglePDF();

        size_t index = m_pdf1.sample(sample.x());

        //warp samples to get vector from triangle using warping function
        Vector3f sampledvector = Warp::squareToUniformTriangle(sample);

        //emitted ray position is found by interpolating the warped directions onto the triangle using the index
        Point3f interpolatedvertex = sampledvector.x() * m_V.col(m_F(0, index)) + sampledvector.y() * m_V.col(m_F(1, index)) + sampledvector.z() * m_V.col(m_F(2, index));
        //Normal3f interpolatednormal = sampledvector.x() * m_N.col(m_F(0, index)) + sampledvector.y() * m_N.col(m_F(1, index)) + sampledvector.z() * m_N.col(m_F(2, index)).normalized();

        eRec.p = interpolatedvertex;

        //computing the normals
        //so if the mesh of triangle gives normal of the point at position p then just interpolate the normal, otherwise find the face normal instead
        if (m_N.size() > 0) {
            eRec.n = sampledvector.x() * m_N.col(m_F(0, index)) + sampledvector.y() * m_N.col(m_F(1, index)) + sampledvector.z() * m_N.col(m_F(2, index)).normalized();
        }
        else {
            Point3f p0 = m_V.col(m_F(0, index));
            Point3f p1 = m_V.col(m_F(1, index));
            Point3f p2 = m_V.col(m_F(2, index));
            Normal3f n = (p1 - p0).cross(p2 - p0).normalized();
            eRec.n = n;
        }
        return eRec.p;
    }

    /// Return a human-readable summary
    std::string toString() const
    {
        return tfm::format(
            "AreaLight[\n"
            "  radiance = %s\n"
            "]", m_radiance.toString());
    }

    EClassType getClassType() const { return EEmitter; }
private:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
