
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>

#include "wprojection.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846264338327
#endif

// Adams TODO list:
// -- Detecting min threshold for kernel extraction prior to bicubic interpolation
//    will be performed by some form of factoring - use of wToMaxSupportRatio and
//    relation to resolution size. eg: 128 (width) and 7 (support)
// -- Some strange pointer bug when using a resolution size < 128

int gridSize = 18000;
double kernelMaxFullSupport = 128.0;
double kernelMinFullSupport = 7.0;

int kernelTextureSize = 16;
int kernelResolutionSize = 128;

double cellsizeRad = 0.000006;
double fieldOfViewDegrees = 0.0;
double maxW = 7000.0;
double wScale = 0.0;
int numWPlanes = 0;
float wToMaxSupportRatio = 0.0;

int main(int argc, char** argv)
{   
    numWPlanes = (int)(maxW * fabs(sin(cellsizeRad * (double)gridSize / 2.0)));
    fieldOfViewDegrees = gridSize * cellsizeRad;
    wScale = pow(numWPlanes - 1, 2.0) / maxW;
    wToMaxSupportRatio = (kernelMaxFullSupport - kernelMinFullSupport) / maxW;
    
    createWProjectionPlanes(kernelResolutionSize, numWPlanes, kernelTextureSize, wScale, fieldOfViewDegrees);
    
    return (EXIT_SUCCESS);
}

int calcWFullSupport(double w, double wToMaxSupportRatio, double minSupport)
{
    // Calculates the full support width of a kernel for w term
    return (int) (fabs(wToMaxSupportRatio * w) + minSupport);
}

void normalizeKernel(DoubleComplex *kernel, int resolution, int support)
{
    // Get sum of magnitudes
    double magnitudeSum;
    int r, c;
    for(r = 0; r < resolution; r++)
        for(c = 0; c < resolution; c++)
            magnitudeSum += complexMagnitude(kernel[r * resolution + c]);
    
    // Normalize weights
    for(r = 0; r < resolution; r++)
        for(c = 0; c < resolution; c++)
            normalizeWeight(&kernel[r * resolution + c], magnitudeSum, resolution, support);
}

DoubleComplex normalizeWeight(DoubleComplex *weight, double mag, int resolution, int support)
{
    weight->real /= mag;
    weight->imaginary /= mag;
    double t2 = (resolution*resolution)/(support*support);
    weight->real *= t2;
    weight->imaginary *= t2;
}

void copyAndTrimKernel(DoubleComplex *dest, DoubleComplex *source, int support)
{
    // Resize and copy over w projection kernel, trim top row and left-most column
    // to ensure center element is the "peak" of the kernel
    DoubleComplex *temp = realloc(dest, (support-1) * (support-1) * sizeof(DoubleComplex));
    if(temp)
    {
        for(int r = 1; r < support; r++)
        {
            for(int c = 1; c < support; c++)
            {   
                int destI = (r-1) * (support-1) + (c-1);
                int srcI = r * support + c;
                dest[destI] = source[srcI];
            }
        }
    }
}

void interpolateKernel(DoubleComplex *source, FloatComplex* dest, int origSupport, int newSupport, int texSupport, int iw, float w, int plane)
{
    int origHalf = (origSupport-1)/2;
    int start = (origHalf/2) * origSupport + (origHalf/2);    
    DoubleComplex *temp = calloc(newSupport * newSupport, sizeof(DoubleComplex));
    
    // Copy critical points of original kernel
    for(int r = 0; r < newSupport; r++)
        for(int c = 0; c < newSupport; c++)
            temp[r * newSupport + c] = source[start + c + (origSupport * r)];
    
    if(iw == plane)
        saveKernelToFile("wproj_%f_trimmed_%d.csv", w, newSupport, temp);
    
    // Perform bicubic interpolation
    InterpolationPoint neighbours[16];
    InterpolationPoint interpolated[4];
    float xShift, yShift, shift2;
    shift2 = getShift(newSupport);
    
    for(int y = 0; y < texSupport; y++)
    {
        yShift = calcShift(y, texSupport);

        for(int x = 0; x < texSupport; x++)
        {
            xShift = calcShift(x, texSupport);
            getBicubicNeighbours(x, y, neighbours, newSupport, texSupport, temp);

            for(int i  = 0, j = -1; i < 4; i++, j++)
            {
                InterpolationPoint newPoint = (InterpolationPoint) {.xShift = xShift, .yShift = (yShift + j * shift2)};
                newPoint = interpolateCubicWeight(neighbours, newPoint, i*4, newSupport, true);
                interpolated[i] = newPoint;
            }

            InterpolationPoint final = (InterpolationPoint) {.xShift = xShift, .yShift = yShift};
            final = interpolateCubicWeight(interpolated, final, 0, newSupport, false);
            int index = y * texSupport + x;
            dest[index] = (FloatComplex) {.real = (float) final.weight.real, .imaginary = (float) final.weight.imaginary};
        }
    }
    
    free(temp);
}

void createWProjectionPlanes(int convolutionSize, int numWPlanes, int textureSupport, double wScale, double fov)
{                
    int convHalf = convolutionSize/2;
    // Flat two dimension array
    DoubleComplex *screen = calloc(convolutionSize * convolutionSize, sizeof(DoubleComplex));
    // Flat cube array (numWPlanes * texFullSupport^2)
    FloatComplex *wTextures = calloc(numWPlanes * textureSupport * textureSupport, sizeof(FloatComplex));
    // Create w screen
    DoubleComplex *shift = calloc(convolutionSize * convolutionSize, sizeof(DoubleComplex));
    // Single dimension spheroidal
    double *spheroidal = calloc(convolutionSize, sizeof(double));
    
    //numWPlanes = 1;
    int plane = numWPlanes-1;
    for(int iw = 0; iw < numWPlanes; iw++)
    {        
        int kernelResolution = convolutionSize;
        // Calculate w term and w specific support size
        double w = iw * iw / wScale;
        int wFullSupport = calcWFullSupport(w, wToMaxSupportRatio, kernelMinFullSupport);
        // Calculate Prolate Spheroidal
        createScaledSpheroidal(spheroidal, wFullSupport, convHalf);
        
        // Create Phase Screen
        createPhaseScreen(kernelResolution, screen, spheroidal, w, fov, wFullSupport);
        
        if(iw == plane)
            saveKernelToFile("wproj_%f_before_fft_%d.csv", w, kernelResolution, screen);
        
        // Perform shift and inverse FFT of Phase Screen
        fft2dShift(kernelResolution, screen, shift);
        inverseFFT2dVectorRadixTransform(kernelResolution, shift, screen);
        fft2dShift(kernelResolution, screen, shift);
        
        if(iw == plane)
            saveKernelToFile("wproj_%f_after_fft_%d.csv", w, kernelResolution, shift);
        
//        // Trim to ensure kernel is odd
//        copyAndTrimKernel(screen, shift, kernelResolution);
//        kernelResolution--;
        // Normalize the kernel
        normalizeKernel(shift, kernelResolution, wFullSupport);
        
        if(iw == plane)
            saveKernelToFile("wproj_%f_normalized_%d.csv", w, kernelResolution, shift);
        
        int support = (kernelResolution - 1) / 2; // half + 1 + half
        FloatComplex *interpolated = calloc(textureSupport * textureSupport, sizeof(FloatComplex));
        interpolateKernel(shift, interpolated, kernelResolution, support, textureSupport, iw, w, plane);
        
        if(iw == plane)
            saveKernelToFile("wproj_%f_interpolated_%d.csv", w, textureSupport, interpolated);
        
        free(interpolated);
        
//        // Locate the support size of kernel stepping in from edge (evaluates magnitude of complex weight)
//        double minThreshold = 0.001;
//        int convHalf = kernelResolution/2;
//        int centerIndex = (kernelResolution*kernelResolution/2) + convHalf;
//        int trialIndex = centerIndex - convHalf;
//        int interpStart;
//        for(interpStart = trialIndex; interpStart < centerIndex; interpStart++)
//        {
//            double magnitude = sqrt(screen[interpStart].real * screen[interpStart].real + screen[interpStart].imaginary * screen[interpStart].imaginary);
//            if(magnitude > minThreshold)
//                break;
//        }
//        
//        int interpHalf = centerIndex - interpStart;
//        int interpWidth = interpHalf * 2 + 1;
//        interpStart -= (kernelResolution * (interpWidth/2));
//        
//        // Allocate memory for sub array and copy
//        DoubleComplex *subArray = calloc(interpWidth * interpWidth, sizeof(DoubleComplex));
//        for(int r = 0; r < interpWidth; r++)
//        {
//            for(int c = 0; c < interpWidth; c++)
//            {
//                int interpIndex = interpStart + c + (kernelResolution * r);
//                subArray[r * interpWidth + c] = screen[interpIndex];
////                if(iw == 0)
////                    printf("%d ", interpIndex);
//            }   
////            if(iw == 0)
////                printf("\n");
//        }
//        
//        InterpolationPoint neighbours[16];
//        InterpolationPoint interpolated[4];
//        
//        float xShift, yShift, shift2;
//        shift2 = getShift(kernelResolution);
//        // Interpolate projection screen into texture
//        for(int y = 0; y < textureSupport; y++)
//        {
//            yShift = calcShift(y, textureSupport);
//            
//            for(int x = 0; x < textureSupport; x++)
//            {
//                xShift = calcShift(x, textureSupport);
//                getBicubicNeighbours(x, y, neighbours, interpWidth, textureSupport, subArray);
//                
//                for(int i  = 0, j = -1; i < 4; i++, j++)
//                {
//                    InterpolationPoint newPoint = (InterpolationPoint) {.xShift = xShift, .yShift = (yShift + j * shift2)};
//                    newPoint = interpolateCubicWeight(neighbours, newPoint, i*4, interpWidth, true);
//                    interpolated[i] = newPoint;
//                }
//                
//                InterpolationPoint final = (InterpolationPoint) {.xShift = xShift, .yShift = yShift};
//                final = interpolateCubicWeight(interpolated, final, 0, interpWidth, false);
//                int index = wTextureIndex(x, y, iw, textureSupport);
//                wTextures[index] = (FloatComplex) {.real = (float) final.weight.real, .imaginary = (float) final.weight.imaginary};
//                if(iw == 0)
//                    printf("%f ", wTextures[index].real);
//            }
//            
//            if(iw == 0)
//                printf("\n");
//        }
//        
//        free(subArray);
        
        // This assumes the reallocate will be successful, consider refactoring
        //realloc(screen, (convolutionSize) * (convolutionSize) * sizeof(DoubleComplex));
        memset(screen, 0, convolutionSize * convolutionSize * sizeof(DoubleComplex));
        memset(shift, 0, convolutionSize * convolutionSize * sizeof(DoubleComplex));
    }
    
    free(spheroidal);
    free(screen);
    free(shift);
    
   
    free(wTextures);
}

void saveKernelToFile(char* filename, float w, int support, DoubleComplex* data)
{
    char *buffer[100];
    sprintf(buffer, filename, w, support);
    FILE *file = fopen(buffer, "w");
    for(int r = 0; r < support; r++)
    {
        for(int c = 0; c < support; c++)
        {
            fprintf(file, "%f, ", data[r * support + c].real);
        }
        fprintf(file, "\n");
    }
    fclose(file);
    printf("FILE SAVED\n");
}

// Col, Row, Depth, Row Dimension
int wTextureIndex(int x, int y, int z, int n)
{
    return ((z)*n*n) + ((y)*n) + (x);
}

void createScaledSpheroidal(double *spheroidal, int wFullSupport, int convHalf)
{
    int wHalfSupport = wFullSupport/2;
    double *nu = calloc(wFullSupport, sizeof(double));
    double *tempSpheroidal = calloc(wFullSupport, sizeof(double));
    // Calculate steps
    for(int i = 0; i < wFullSupport; i++)
        nu[i] = fabs(calcShift(i, wFullSupport));
    // Calculate curve from steps
    calcSpheroidalCurve(nu, tempSpheroidal, wFullSupport);
    // Bind weights to middle
    for(int i = convHalf-wHalfSupport; i <= convHalf+wHalfSupport; i++)
        spheroidal[i] = tempSpheroidal[i-(convHalf-wHalfSupport)];
    free(tempSpheroidal);
    free(nu);
}

void createPhaseScreen(int convSize, DoubleComplex *screen, double* spheroidal, double w, double fieldOfView, int scalarSupport)
{        
    int convHalf = convSize/2;
    int scalarHalf = scalarSupport/2;
    int index = 0;
    double taper, taperY;
    double l, m;
    double lsq, rsq;
    double phase;
    
    int panicCounter = 0;
    
    for(int iy = 0; iy < scalarSupport; iy++)
    {
        l = (((double) iy-(scalarHalf+0.5)) / (double) scalarSupport) * fieldOfView;
        lsq = l*l;
        taperY = spheroidal[iy+(convHalf-scalarHalf)];
        phase = 0.0;
        
        for(int ix = 0; ix < scalarSupport; ix++)
        {
            m = (((double) ix-(scalarHalf+0.5)) / (double) scalarSupport) * fieldOfView;
            rsq = lsq+(m*m);
            taper = taperY * spheroidal[ix+(convHalf-scalarHalf)];
            index = (iy+(convHalf-scalarHalf)) * convSize + (ix+(convHalf-scalarHalf));
            
            if(rsq < 1.0)
            {
                phase = w * (1.0 - sqrt(1.0 - rsq));
                screen[index] = complexConjugateExp(phase);
            }
            
            if(rsq == 0.0)
                screen[index] = (DoubleComplex) {.real = 1.0, .imaginary = 0.0};
            
            // Bug fix: what happens when W = 0?
                
            if(rsq >= 1.0)
                panicCounter++;
                
            // printf("%f ", taper);
                
            screen[index].real *= taper;
            screen[index].imaginary *= taper;
        }
        // printf("\n");
    }
    if(panicCounter > 0)
        printf("You have panicked a total of %d times, god damn!\n", panicCounter);
}

void inverseFFT2dVectorRadixTransform(int numChannels, DoubleComplex *input, DoubleComplex *output)
{   
    // Calculate bit reversed indices
    int* bitReversedIndices = malloc(sizeof(int) * numChannels);
    calcBitReversedIndices(numChannels, bitReversedIndices);
    
    // Copy data to result for processing
    for(int r = 0; r < numChannels; r++)
        for(int c = 0; c < numChannels; c++)
            output[r * numChannels + c] = input[bitReversedIndices[r] * numChannels + bitReversedIndices[c]];
    free(bitReversedIndices);
    
    // Use butterfly operations on result to find the DFT of original data
    for(int m = 2; m <= numChannels; m *= 2)
    {
        DoubleComplex omegaM = (DoubleComplex) {.real = cos(M_PI * 2.0 / m), .imaginary = sin(M_PI * 2.0 / m)};
        
        for(int k = 0; k < numChannels; k += m)
        {
            for(int l = 0; l < numChannels; l += m)
            {
                DoubleComplex x = (DoubleComplex) {.real = 1.0, .imaginary = 0.0};
                
                for(int i = 0; i < m / 2; i++)
                {
                    DoubleComplex y = (DoubleComplex) {.real = 1.0, .imaginary = 0.0};
                    
                    for(int j = 0; j < m / 2; j++)
                    {   
                        // Perform 2D butterfly operation in-place at (k+j, l+j)
                        int in00Index = (k+i) * numChannels + (l+j);
                        DoubleComplex in00 = output[in00Index];
                        int in01Index = (k+i) * numChannels + (l+j+m/2);
                        DoubleComplex in01 = complexMultiply(output[in01Index], y);
                        int in10Index = (k+i+m/2) * numChannels + (l+j);
                        DoubleComplex in10 = complexMultiply(output[in10Index], x);
                        int in11Index = (k+i+m/2) * numChannels + (l+j+m/2);
                        DoubleComplex in11 = complexMultiply(complexMultiply(output[in11Index], x), y);
                        
                        DoubleComplex temp00 = complexAdd(in00, in01);
                        DoubleComplex temp01 = complexSubtract(in00, in01);
                        DoubleComplex temp10 = complexAdd(in10, in11);
                        DoubleComplex temp11 = complexSubtract(in10, in11);
                        
                        output[in00Index] = complexAdd(temp00, temp10);
                        output[in01Index] = complexAdd(temp01, temp11);
                        output[in10Index] = complexSubtract(temp00, temp10);
                        output[in11Index] = complexSubtract(temp01, temp11);
                        y = complexMultiply(y, omegaM);
                    }
                    x = complexMultiply(x, omegaM);
                }
            }
        }
    }
    
    for(int i = 0; i < numChannels; i++)
        for(int j = 0; j < numChannels; j++)
        {
            output[i * numChannels + j].real /= (numChannels * numChannels);
            output[i * numChannels + j].imaginary /= (numChannels * numChannels);
        }
}

void calcBitReversedIndices(int n, int* indices)
{   
    for(int i = 0; i < n; i++)
    {
        // Calculate index r to which i will be moved
        unsigned int iPrime = i;
        int r = 0;
        for(int j = 1; j < n; j*=2)
        {
            int b = iPrime & 1;
            r = (r << 1) + b;
            iPrime = (iPrime >> 1);
        }
        indices[i] = r;
    }
}

void calcSpheroidalCurve(double *nu, double *curve, int width)
{   
    double p[2][5] = {{8.203343e-2, -3.644705e-1, 6.278660e-1, -5.335581e-1, 2.312756e-1},
                     {4.028559e-3, -3.697768e-2, 1.021332e-1, -1.201436e-1, 6.412774e-2}};
    double q[2][3] = {{1.0000000e0, 8.212018e-1, 2.078043e-1},
                     {1.0000000e0, 9.599102e-1, 2.918724e-1}};
    
    int pNum = 5;
    int qNum = 3;
    
    int *part = calloc(width, sizeof(int));
    double *nuend = calloc(width, sizeof(double));
    double *delnusq = calloc(width, sizeof(double));
    double *top = calloc(width, sizeof(double));
    double *bottom = calloc(width, sizeof(double));
    
    for(int i = 0; i < width; i++)
    {
        if(nu[i] >= 0.0 && nu[i] <= 0.75)
            part[i] = 0;
        if(nu[i] > 0.75 && nu[i] < 1.0)
            part[i] = 1;
        
        if(nu[i] >= 0.0 && nu[i] <= 0.75)
            nuend[i] = 0.75;
        if(nu[i] > 0.75 && nu[i] < 1.0)
            nuend[i] = 1.0;
        
        delnusq[i] = (nu[i] * nu[i]) - (nuend[i] * nuend[i]);
    }
    
    for(int i = 0; i < width; i++)
    {
        top[i] = p[part[i]][0];
        bottom[i] = q[part[i]][0];
    }
    
    for(int i = 1; i < pNum; i++)
        for(int y = 0; y < width; y++)
            top[y] += (p[part[y]][i] * pow(delnusq[y], i)); 
    
    for(int i = 1; i < qNum; i++)
        for(int y = 0; y < width; y++)
            bottom[y] += (q[part[y]][i] * pow(delnusq[y], i));
    
    for(int i = 0; i < width; i++)
    {   
        if(bottom[i] > 0.0)
            curve[i] = top[i] / bottom[i];
        if(fabs(nu[i]) > 1.0)
            curve[i] = 0.0;
        
    }
    
    free(bottom);
    free(top);
    free(delnusq);
    free(nuend);
    free(part);
}

DoubleComplex complexAdd(DoubleComplex x, DoubleComplex y)
{
    DoubleComplex z;
    z.real = x.real + y.real;
    z.imaginary = x.imaginary + y.imaginary;
    return z;
}

DoubleComplex complexSubtract(DoubleComplex x, DoubleComplex y)
{
    DoubleComplex z;
    z.real = x.real - y.real;
    z.imaginary = x.imaginary - y.imaginary;
    return z;
}

DoubleComplex complexMultiply(DoubleComplex x, DoubleComplex y)
{
    DoubleComplex z;
    z.real = x.real*y.real - x.imaginary*y.imaginary;
    z.imaginary = x.imaginary*y.real + x.real*y.imaginary;
    return z;
}

double complexMagnitude(DoubleComplex x)
{
    return sqrt(x.real * x.real + x.imaginary * x.imaginary);
}

DoubleComplex complexDivide(DoubleComplex x, DoubleComplex y) 
{
    DoubleComplex z;
    double denominator = (y.real * y.real) + (y.imaginary * y.imaginary);  
    z.real = (x.real*y.real + x.imaginary*y.imaginary) / denominator;
    z.imaginary = (x.imaginary*y.real - x.real*y.imaginary) / denominator;
    return z;    
}

DoubleComplex complexConjugateExp(double ph)
{
    return (DoubleComplex) {.real = cos((double)(2.0*M_PI*ph)), .imaginary = -sin((double)(2.0*M_PI*ph))};
}

void fft2dShift(int n, DoubleComplex *input, DoubleComplex *shifted)
{
    // Can be refactored to save memory
    int r = 0, c = 0;
    for(int i = -n/2; i < n/2; i++)
    {
        for(int j = -n/2; j < n/2; j++)
        {
            if(i >= 0 && j >= 0)
                shifted[r * n + c] = input[i * n +j];
            else if(i < 0 && j >=0)
                shifted[r * n + c] = input[(i+n)*n+j];
            else if(i >= 0 && j < 0)
                shifted[r * n + c] = input[i*n+j+n];
            else
                shifted[r * n + c] = input[(i+n)*n+j+n];
            
            c++;
        }
        r++;
        c = 0;
    }
}

InterpolationPoint interpolateCubicWeight(InterpolationPoint *points, InterpolationPoint newPoint, int start, int width, bool horizontal)
{      
    // Convert to double?
    double shiftCubed = pow(getShift(width), 3);

    DoubleComplex p0 = (DoubleComplex) {.real = (horizontal) ? points[start+0].xShift : points[start+0].yShift, .imaginary = 0.0};
    DoubleComplex p1 = (DoubleComplex) {.real = (horizontal) ? points[start+1].xShift : points[start+1].yShift, .imaginary = 0.0};
    DoubleComplex p2 = (DoubleComplex) {.real = (horizontal) ? points[start+2].xShift : points[start+2].yShift, .imaginary = 0.0};
    DoubleComplex p3 = (DoubleComplex) {.real = (horizontal) ? points[start+3].xShift : points[start+3].yShift, .imaginary = 0.0};
    DoubleComplex interpShift = (DoubleComplex) {.real = (horizontal) ? newPoint.xShift : newPoint.yShift, .imaginary = 0.0};
    
    DoubleComplex w0 = (DoubleComplex) {.real = -(points[start+0].weight.real) / (6.0 * shiftCubed), 
            .imaginary = -(points[start+0].weight.imaginary) / (6.0 * shiftCubed)};
    DoubleComplex w1 = (DoubleComplex) {.real = points[start+1].weight.real / (2.0 * shiftCubed),
            .imaginary = points[start+1].weight.imaginary / (2.0 * shiftCubed)};
    DoubleComplex w2 = (DoubleComplex) {.real = -points[start+2].weight.real / (2.0 * shiftCubed), 
            .imaginary = -points[start+2].weight.imaginary / (2.0 * shiftCubed)};
    DoubleComplex w3 = (DoubleComplex) {.real = points[start+3].weight.real / (6.0 * shiftCubed), 
            .imaginary = points[start+3].weight.imaginary / (6.0 * shiftCubed)}; 
    
    // Refactor for complex multiplication and subtraction
    DoubleComplex t0 = complexMultiply(complexMultiply(complexMultiply(w0, complexSubtract(interpShift, p1)), complexSubtract(interpShift, p2)), 
            complexSubtract(interpShift, p3));
    DoubleComplex t1 = complexMultiply(complexMultiply(complexMultiply(w1, complexSubtract(interpShift, p0)), complexSubtract(interpShift, p2)), 
            complexSubtract(interpShift, p3));
    DoubleComplex t2 = complexMultiply(complexMultiply(complexMultiply(w2, complexSubtract(interpShift, p0)), complexSubtract(interpShift, p1)),
            complexSubtract(interpShift, p3));
    DoubleComplex t3 = complexMultiply(complexMultiply(complexMultiply(w3, complexSubtract(interpShift, p0)), complexSubtract(interpShift, p1)), 
            complexSubtract(interpShift, p2));
    // Refactor for complex addition
    newPoint.weight = complexAdd(complexAdd(complexAdd(t0, t1), t2), t3);
    return newPoint;
}

void getBicubicNeighbours(int x, int y, InterpolationPoint *neighbours, int kernelFullSupport, int interpFullSupport, DoubleComplex* matrix)
{
    // Transform x, y into scaled shift
    float shiftX = calcShift(x, interpFullSupport);
    float shiftY = calcShift(y, interpFullSupport);
    // Get x, y from scaled shift 
    int scaledPosX = calcPosition(shiftX, kernelFullSupport);
    int scaledPosY = calcPosition(shiftY, kernelFullSupport);
    // Get 16 neighbours
    for(int r = scaledPosY - 1, i = 0; r < scaledPosY + 3; r++)
    {
        for(int c = scaledPosX - 1; c < scaledPosX + 3; c++)
        {
            InterpolationPoint n = (InterpolationPoint) {.xShift = calcShift(c, kernelFullSupport), .yShift = calcShift(r, kernelFullSupport)};
            
            if(c < 0 || c > kernelFullSupport || r < 0 || r > kernelFullSupport)
                n.weight = (DoubleComplex) {.real = 0.0, .imaginary = 0.0};
            else
                n.weight = matrix[r * kernelFullSupport + c];
            
            neighbours[i++] = n;
        }
    }
}

float calcShift(int index, int width)
{
    return getStartShift(width) + index * getShift(width);
}

double getShift(double width)
{
    return 2.0/width;
}

float getStartShift(float width)
{
    return -1.0 + (1.0 / width);
}

int calcPosition(float x, int scalerWidth)
{
    return (int) floor(((x+1.0f)/2.0f) * (scalerWidth));
}

