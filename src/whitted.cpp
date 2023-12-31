#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator
{
public:
    WhittedIntegrator(const PropertyList& props)
    {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {

        
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(1.0f);

        Color3f Le(0.0f);
        if (its.mesh->isEmitter())
        {
            EmitterQueryRecord eRec(ray.o, its.p, its.shFrame.n);
            Le = its.mesh->getEmitter()->eval(eRec);

        }

        //Intersection point on scene where we want to calculate radiance
        Point3f x = its.p;
        const Mesh* intersectedmesh = its.mesh; //mesh of the scene
        const BSDF* bsdf = intersectedmesh->getBSDF(); //bsdf of the scene

        //if it is a diffuse material
        if (bsdf->isDiffuse())
        {
            //get all the meshes of the emitter in a vector
            std::vector<Mesh*> allmeshes = scene->getEmittermeshes();

            //sample randomly to get one of the meshes to be emitted from
            int numberofmeshes = allmeshes.size();
            float randomemittersample = (sampler->next1D()) * numberofmeshes;
            int emittertosamplefrom = (int)(randomemittersample);

            //get the mesh of this emittertosamplefrom
            Mesh* emittertoemit = allmeshes[emittertosamplefrom];

            //get the emitter instance of the emitter we have chosen
            Emitter* emittermesh = emittertoemit->getEmitter();

            //Make a structure of emitter having origin at x  
            EmitterQueryRecord EQR(x);

            //根据光源采样
            Point3f y = emittermesh->sample(emittertoemit, EQR, sampler);
            Color3f emitterradiance = emittermesh->eval(EQR);
            float pdf = emittermesh->pdf(emittertoemit,EQR);

            //checking visibility of y and x
            Vector3f dir = y - x;
            Vector3f oppdir = x - y;
            Ray3f shadowRay = Ray3f(x, dir);
            Point3f lightnormal = EQR.n;

            //checks for emitter point sampled normal to always be in downward direction
            float check = lightnormal.dot(oppdir.normalized());
            //if it is in upward direction; then ignore the ray
            if (check < 0.0f)
                return Color3f(0.0f);
            else
            {
                Intersection its1;
                int v = 0;
                if (scene->rayIntersect(shadowRay, its1))
                {
                    if (its1.mesh == emittertoemit)
                        v = 1;
                }


                //finding costheta value; that is normal at x and ray to emitter
                float costheta = std::abs(its.shFrame.cosTheta(its.shFrame.toLocal(dir).normalized()));
                //finding cosalpha value; that is normal at y and ray to point x
                float cosalpha = std::abs(lightnormal.dot(oppdir.normalized()));
                //Calculating the bsdf

                //getting the origin point of ray from camera from w_o
                Vector3f w_o = its.toLocal((ray.o - x).normalized());
                BSDFQueryRecord bsdfquery(its.shFrame.toLocal(dir), w_o, EMeasure::ESolidAngle);
                //get bsdf for the scene at point x
                Color3f brdfofmesh = bsdf->eval(bsdfquery);

                float r = (brdfofmesh.r() * costheta * cosalpha * emitterradiance.r() * v) / (pdf * dir.norm() * dir.norm());
                float g = (brdfofmesh.g() * costheta * cosalpha * emitterradiance.g() * v) / (pdf * dir.norm() * dir.norm());
                float b = (brdfofmesh.b() * costheta * cosalpha * emitterradiance.b() * v) / (pdf * dir.norm() * dir.norm());
                //cout << (pdf * dir.norm() * dir.norm()) <<" " << r << " " << g << " " << b << endl;
                //exit(0);
                return Color3f(r, g, b) + Le;
            }
        }
        //if the material bsdf is not diffuse
        else
        {
            Vector3f w_o = its.toLocal((ray.o - x).normalized());
            BSDFQueryRecord bsdfquery(w_o);
            //get bsdf for scene at point x

            Color3f brdfofmesh = bsdf->sample(bsdfquery, sampler->next2D());

            float addsampler = sampler->next1D();
            Ray3f newray = Ray3f(x, its.toWorld(bsdfquery.wo));
            if (addsampler < 0.95f)
            {
                Color3f rad;
                rad = (1.0f / 0.95f) * (brdfofmesh)*Li(scene, sampler, newray);
                return rad;
            }
            else
            {
                return Color3f(0.0f) + Le;
            }
        }
    }

    std::string toString() const
    {
        return tfm::format(
            "WhittedIntegrator[]"
        );
    }
protected:
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END
