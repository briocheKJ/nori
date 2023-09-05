#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/common.h>
#include <Python.h>

NORI_NAMESPACE_BEGIN



class NeuralBRDF : public BSDF {
public:
    NeuralBRDF(const PropertyList& propList) {
        Py_SetPythonHome(L"C:\\Users\\brioche\\anaconda3");

        // 初始化Python解释器
        Py_Initialize();

        // 导入模块
        PyObject* pName = PyUnicode_DecodeFSDefault("hello");
        pModule = PyImport_Import(pName); 
        Py_XDECREF(pName);

        if (pModule != NULL) {
            // 获取函数引用
            PyObject* pFunc1 = PyObject_GetAttrString(pModule, "load_model");

            if (PyCallable_Check(pFunc1)) {
                // 调用函数
                PyObject* pValue = PyObject_CallObject(pFunc1, NULL);

                // TODO: 检查和使用 pValue (如果需要)

                Py_XDECREF(pValue);
            }
            else {
                PyErr_Print();
            }
            Py_XDECREF(pFunc);
            Py_XDECREF(pModule);
        }
        else {
            PyErr_Print();
        }

        if (pModule != NULL) {
            pFunc = PyObject_GetAttrString(pModule, "compute_rgb");
            if (!PyCallable_Check(pFunc)) {
                PyErr_Print();
                Py_XDECREF(pFunc);
                pFunc = NULL;
            }
        }
        else {
            PyErr_Print();
        }
        PyEval_InitThreads();
        PyThreadState* mainThreadState = PyEval_SaveThread();

        /*PyObject* pFunc;
        PyObject* pArgs, * pValue;
        pName = PyUnicode_DecodeFSDefault("hello");  // 指定Python模块名
        pModule = PyImport_Import(pName);
        Py_XDECREF(pName);

        if (pModule != NULL) {
            pFunc = PyObject_GetAttrString(pModule, "load_model");  // 指定函数名
            if (PyCallable_Check(pFunc)) {
                pArgs = PyTuple_Pack(1, PyLong_FromLong(123));  // 传递一个参数：整数123
                pValue = PyObject_CallObject(pFunc, pArgs);
                Py_XDECREF(pArgs);
                if (pValue != NULL) {
                    printf("Returned value: %ld\n", PyLong_AsLong(pValue));
                    Py_XDECREF(pValue);
                }
            }
            else {
                PyErr_Print();
            }
            Py_XDECREF(pFunc);
            Py_XDECREF(pModule);
        }
        else {
            PyErr_Print();
        }*/
        //cout << 1 << endl;
    }


    // Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord& bRec) const {

        Vector3f wi = bRec.wi.normalized();
        Vector3f wo = bRec.wo.normalized();
        // 导入模块
        float result1, result2, result3;
        PyGILState_STATE state = PyGILState_Ensure();
        PyObject* pArgs = PyTuple_Pack(6, PyFloat_FromDouble(wi.x()), PyFloat_FromDouble(wi.y()),
            PyFloat_FromDouble(wi.z()), PyFloat_FromDouble(wo.x()),
            PyFloat_FromDouble(wo.y()), PyFloat_FromDouble(wo.z()));
       
        PyObject* pValue = PyObject_CallObject(pFunc, pArgs);
        Py_XDECREF(pArgs);

        if (pValue) {
            result1 = PyFloat_AsDouble(PyTuple_GetItem(pValue, 0));
            result2 = PyFloat_AsDouble(PyTuple_GetItem(pValue, 1));
            result3 = PyFloat_AsDouble(PyTuple_GetItem(pValue, 2));
            Py_XDECREF(pValue);
        }
        else {
            PyErr_Print();
        } 
        PyGILState_Release(state);
        //cout << result1 << " " << result2 << " " << result3 << endl;
        //exit(0);
        return Color3f(result1,result2,result3);
    }


    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord& bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;


        /* Importance sampling density wrt. solid angles:
           cos(theta) / pi.

           Note that the directions in 'bRec' are in local coordinates,
           so Frame::cosTheta() actually just returns the 'z' component.
        */
        return INV_PI * Frame::cosTheta(bRec.wo);
    }


    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord& bRec, const Point2f& sample) const {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;

        /* Warp a uniformly distributed sample on [0,1]^2
           to a direction on a cosine-weighted hemisphere */
        bRec.wo = Warp::squareToCosineHemisphere(sample);

        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;

        /* eval() / pdf() * cos(theta) = albedo. There
           is no need to call these functions. */
        return Color3f(1.0);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "NeuralBRDF[]"
        );
    }
private:
    PyObject* pModule;
    PyObject* pFunc = NULL;
};

NORI_REGISTER_CLASS(NeuralBRDF, "neuralBRDF");
NORI_NAMESPACE_END

