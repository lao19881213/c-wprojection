
#ifndef WPROJECTION_H
#define WPROJECTION_H

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct FloatComplex {
        float real;
        float imaginary;
    } FloatComplex;

    typedef struct DoubleComplex {
        double real;
        double imaginary;
    } DoubleComplex;

    typedef struct InterpolationPoint {
        float xShift;
        float yShift;
        DoubleComplex weight;
    } InterpolationPoint;

    double getShift(double width);
    float getStartShift(float width);
    float calcShift(int index, int width, float start);
    int calcPosition(float x, int scalerWidth);
    InterpolationPoint interpolateCubicWeight(InterpolationPoint *points, InterpolationPoint newPoint, int start, int width, bool horizontal);
    void getBicubicNeighbours(int x, int y, InterpolationPoint *neighbours, int kernelFullSupport, int interpFullSupport, DoubleComplex* matrix);

    void createScaledSpheroidal(double *spheroidal, int wFullSupport, int convHalf);
    int calcWFullSupport(double w, double wToMaxSupportRatio, double minSupport);
    void createWProjectionPlanes(int convolutionSize, int numWPlanes, int textureSupport, double wScale, double fov);
    void createPhaseScreen(int convSize, DoubleComplex *wScreen, double* spheroidal, double w, double fieldOfView, int scalarSupport);
    void calcSpheroidalCurve(double *nu, double *curve, int width);
    void inverseFFT2dVectorRadixTransform(int numChannels, DoubleComplex *input, DoubleComplex *output);
    void calcBitReversedIndices(int n, int* indices);
    void fft2dShift(int n, DoubleComplex *input, DoubleComplex *shifted);
    float calcSpheroidalShift(int index, int width);
    float calcInterpolateShift(int index, int width, float start);

    void saveKernelToFile(char* filename, float w, int support, DoubleComplex* data);
    void copyAndTrimKernel(DoubleComplex *dest, DoubleComplex *source, int support);
    void normalizeKernel(DoubleComplex *kernel, int resolution, int support);
    DoubleComplex normalizeWeight(DoubleComplex weight, double mag, int resolution, int support);
    void interpolateKernel(DoubleComplex *source, DoubleComplex* dest, int origSupport, int texSupport);

    DoubleComplex complexAdd(DoubleComplex x, DoubleComplex y);
    DoubleComplex complexSubtract(DoubleComplex x, DoubleComplex y);
    DoubleComplex complexMultiply(DoubleComplex x, DoubleComplex y);
    DoubleComplex complexConjugateExp(double ph);
    double complexMagnitude(DoubleComplex x);

#ifdef __cplusplus
}
#endif

#endif /* WPROJECTION_H */
