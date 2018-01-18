
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

int gridSize = 18000;
double kernelMaxFullSupport = 128.0;
double kernelMinFullSupport = 7.0;

int kernelTextureSize = 256;
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
    
    printf("Num Planes: %f\n",wScale);
    
    createWProjectionPlanes(kernelResolutionSize, numWPlanes, kernelTextureSize, wScale, fieldOfViewDegrees, wToMaxSupportRatio, kernelMinFullSupport);
    
    return (EXIT_SUCCESS);
}

double calcWFullSupport(double w, double wToMaxSupportRatio, double minSupport)
{
    // Calculates the full support width of a kernel for w term
    return fabs(wToMaxSupportRatio * w) + minSupport;
}

void createWProjectionPlanes(int convolutionSize, int numWPlanes, int textureSupport, double wScale, double fov,
        double wToMaxSupportRatio, double kernelMinSupport)
{                
    double *nu = calloc(convolutionSize, sizeof(double));
    double *spheroidal = calloc(convolutionSize, sizeof(double));
    // Calculate steps
    for(int i = 0; i < convolutionSize; i++)
        nu[i] = fabs(calcShift(i, convolutionSize));
        
    // Calculate curve from steps
    calcSpheroidalCurve(nu, spheroidal, convolutionSize);
    
    double sphrMax = -DBL_MAX;
    for(int r = 0; r < convolutionSize; r++)
        for(int c = 0; c < convolutionSize; c++)
            if(spheroidal[r] * spheroidal[c] > sphrMax)
                sphrMax = spheroidal[r] * spheroidal[c];
    free(nu);   
    
    // Flat two dimension array
    DoubleComplex *screen = calloc(convolutionSize * convolutionSize, sizeof(DoubleComplex));
    // Flat cube array (numWPlanes * texFullSupport^2)
    FloatComplex *wTextures = calloc(numWPlanes * textureSupport * textureSupport, sizeof(FloatComplex));
    // Create w screen
    DoubleComplex *shift = calloc(convolutionSize * convolutionSize, sizeof(DoubleComplex));
    
    // numWPlanes = 1;
    for(int iw = 0; iw < numWPlanes; iw++)
    {        
        double w = iw * iw / wScale;
        double fresnel = w * ((0.5 * fov)*(0.5 * fov));
        printf("CreateWTermLike: For w = %f, field of view = %f, fresnel number = %f\n", w, fov, fresnel);
        double fovScale = 1.0;//(double) convolutionSize / calcWFullSupport(w, wToMaxSupportRatio, kernelMinSupport);
        printf("FOV Scale: %f\n", fovScale);
        createPhaseScreen(convolutionSize, screen, spheroidal, w, fov, sphrMax, fovScale);
        
        // Shift FFT for transformation
        fft2dShift(convolutionSize, screen, shift);
        inverseFFT2dVectorRadixTransform(convolutionSize, shift, screen);
        fft2dShift(convolutionSize, screen, shift);
        
                
        if(iw == 0)
        {
            // Print to file
            char *buffer[100];
            sprintf(buffer, "wproj_screen_blah_%d.csv", convolutionSize);
            FILE *file = fopen(buffer, "w");
            for(int r = 0; r < convolutionSize; r++)
            {
                for(int c = 0; c < convolutionSize; c++)
                {
                    fprintf(file, "%f, ", shift[r * convolutionSize + c].real);
                }
                fprintf(file, "\n");
            }
            fclose(file);
            printf("FILE SAVED\n");
        }
        
        // Reallocate wKernel for mirroring missing row/cols
        DoubleComplex *temp = realloc(screen, (convolutionSize+1) * (convolutionSize+1) * sizeof(DoubleComplex));
        if(temp)
            // does temp become NULL?
            screen = temp;
        else
            printf("ERR: Unable to reallocate wKernel\n");
        
        // Append to list of kernels
        for(int r = 0; r < convolutionSize + 1; r++)
        {
            for(int c = 0; c < convolutionSize + 1; c++)
            {   
                // Copy existing weights
                if(r < convolutionSize && c < convolutionSize)
                    screen[r * convolutionSize + c] = shift[r * convolutionSize + c];
                // Mirror top row onto bottom row
                if(r == convolutionSize && c < convolutionSize)
                    screen[r * convolutionSize + c] = screen[0 * convolutionSize + c];
                // Mirror far left col weight onto far right col
                if(c == convolutionSize)
                    screen[r * convolutionSize + c] = screen[r * convolutionSize + 0];
            }
        }
        
        
        
        
        // Locate the support size of kernel stepping in from edge (evaluates magnitude of complex weight)
        double minThreshold = 0.001;
        int convHalf = convolutionSize/2;
        int centerIndex = (convolutionSize*convolutionSize/2) + convHalf;
        printf("Center: %d\n", centerIndex);
        int trialIndex = centerIndex - convHalf;
        int interpStart;
        for(interpStart = trialIndex; interpStart < centerIndex; interpStart++)
        {
            double magnitude = screen[interpStart].real;//sqrt(screen[interpStart].real * screen[interpStart].real + screen[interpStart].imaginary * screen[interpStart].imaginary);
            //printf("%d: Val: %f\n", interpStart, magnitude);
            if(magnitude > minThreshold)
                break;
        }
        
        int interpHalf = centerIndex - interpStart;
        int interpWidth = interpHalf * 2 + 1;
        interpStart -= (convolutionSize * (interpWidth/2));
        
        // Allocate memory for sub array and copy
        DoubleComplex *subArray = calloc(interpWidth * interpWidth, sizeof(DoubleComplex));
        for(int r = 0; r < interpWidth; r++)
        {
            for(int c = 0; c < interpWidth; c++)
            {
                int interpIndex = interpStart + c + (convolutionSize * r);
                subArray[r * interpWidth + c] = screen[interpIndex];
//                if(iw == 0)
//                    printf("%d ", interpIndex);
            }   
//            if(iw == 0)
//                printf("\n");
        }
        
        InterpolationPoint neighbours[16];
        InterpolationPoint interpolated[4];
        
        float xShift, yShift, shift2;
        shift2 = getShift(convolutionSize);
        // Interpolate projection screen into texture
        for(int y = 0; y < textureSupport; y++)
        {
            yShift = calcShift(y, textureSupport);
            
            for(int x = 0; x < textureSupport; x++)
            {
                xShift = calcShift(x, textureSupport);
                getBicubicNeighbours(x, y, neighbours, interpWidth, textureSupport, subArray);
                
                for(int i  = 0, j = -1; i < 4; i++, j++)
                {
                    InterpolationPoint newPoint = (InterpolationPoint) {.xShift = xShift, .yShift = (yShift + j * shift2)};
                    newPoint = interpolateCubicWeight(neighbours, newPoint, i*4, interpWidth, true);
                    interpolated[i] = newPoint;
                }
                
                InterpolationPoint final = (InterpolationPoint) {.xShift = xShift, .yShift = yShift};
                final = interpolateCubicWeight(interpolated, final, 0, interpWidth, false);
                int index = wTextureIndex(x, y, iw, textureSupport);
                wTextures[index] = (FloatComplex) {.real = (float) final.weight.real, .imaginary = (float) final.weight.imaginary};
//                if(iw == 0)
//                    printf("%f ", wTextures[index].real);
            }
            
//            if(iw == 0)
//                printf("\n");
        }
        
        free(subArray);
        
        
//        // Trim kernels to texture dimensions
//        int kernelCenter = (convolutionSize / 2);
//        for(int r = 0; r < textureSupport; r++)
//        {
//            int kernelRowIndex = (kernelCenter + r) - (textureSupport/2);
//
//            for(int c = 0; c < textureSupport; c++)
//            {
//                int kernelColIndex = (kernelCenter + c) - (textureSupport/2);
//                int index = wTextureIndex(c, r, iw, textureSupport);
//                // Normal assignment
//                wTextures[index] = (FloatComplex) {.real = (float) screen[kernelRowIndex * convolutionSize + kernelColIndex].real,
//                    .imaginary = (float) screen[kernelRowIndex * convolutionSize + kernelColIndex].imaginary};
//
////                printf("%4.10f%+4.10f ", wTextures[index].real, wTextures[index].imaginary);
//
//            }
//        }
        
        memset(screen, 0, convolutionSize * convolutionSize * sizeof(DoubleComplex));
        memset(shift, 0, convolutionSize * convolutionSize * sizeof(DoubleComplex));
    }
    
    free(spheroidal);
    free(screen);
    free(shift);
    
//    for(int iw = 0; iw < numWPlanes; iw++)
//    {
//        for(int iy = 0; iy < textureSupport; iy++)
//        {
//            for(int ix = 0; ix < textureSupport; ix++)
//            {
//                int i = wTextureIndex(ix, iy, iw, textureSupport);
//                
//                if(iw == 0)
//                    printf("%.20f, ", wTextures[i].real);
//            }
//            
//            if(iw == 0)
//                printf("\n");
//        }
//    }
    
    free(wTextures);
}

// Col, Row, Depth, Row Dimension
int wTextureIndex(int x, int y, int z, int n)
{
    return ((z)*n*n) + ((y)*n) + (x);
}

void createPhaseScreen(int convSize, DoubleComplex *screen, double* spheroidal, double w, double fieldOfView, double taperMax, float fovScale)
{        
    int convHalf = convSize/2;
    int index = 0;
    double taper, taperY;
    double l, m;
    double lsq, rsq;
    double phase;
    
    int panicCounter = 0;
    
    for(int iy = 0; iy < convSize; iy++)
    {
        l = (((double) iy-(convHalf+0.5)) / (double) convSize) * fieldOfView * fovScale;
        lsq = l*l;
        taperY = spheroidal[iy];
        phase = 0.0;
        
        for(int ix = 0; ix < convSize; ix++)
        {
            m = (((double) ix-(convHalf+0.5)) / (double) convSize) * fieldOfView * fovScale;
            rsq = lsq+(m*m);
            taper = taperMax / (taperY * spheroidal[ix]);
            index = iy * convSize + ix;
            
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
                
            screen[index].real /= taper;
            screen[index].imaginary /= taper;
        }
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

