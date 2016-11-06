// Base code for 2D FFT


// Globals

int N;                          // image dimensions NxN.  N is power of two.
PImage img;                     // image read from a file, padded to NxN
float[] spatialSignal;          // spatial-domain image intensities, size NxN
Complex[] frequencySignal;      // frequency-domain image, size NxN

PFont font;                     // for drawing text
PImage outputImage;             // image currently being rendered, size NxN

boolean showingFT = false;      // true when showing the F.T.

int filterType = 0;             // type of filter to apply to F.T. (0=Gaussian, 1=Butterworth, 2=box)
float filterRadius = 15;        // radius of filter applied to image (in units of pixels)

float logFactor = 400;          // log factor for enhancing spectrum


// "constants"

int winWidth = 1024;            // window width
int winHeight = 1024;           // window height

int minLum = 16;                // min luminance (Y) value storable
int maxLum = 235;               // max luminance (Y) value storable

String imageFilename = "noisy1.jpg"; // current image file


// Init


void settings() {

    // Set up screen
    
    size( winWidth, winHeight, P2D );
    smooth( 2 );
}


void setup()

{
    // Set up text font

    font = createFont( "Arial", 16, true );
    textFont( font );

    // load default image upon startup

    readImage();
}


// Draw to window

void draw() {

    background( color(255,255,255) );

    // Show the image

    image( outputImage, 0, 0 );

    // Show the mouse

    if (showingFT) {
        noFill();
        stroke( 255 );
        ellipse( mouseX, mouseY, filterRadius*2, filterRadius*2 );
    }

    // Status message: "image filename | F.T. enhancement | F.T. filter"

    String msg = imageFilename;

    if (showingFT) {

        msg = msg + "  |  enhancement " + nf( logFactor, 0, 1 ) + "  |  ";

        switch (filterType) {
        case 0: msg = msg + "Gaussian filter"; break;
        case 1: msg = msg + "Butterworth filter"; break;
        case 2: msg = msg + "box filter"; break;
        }

        msg = msg + ", radius " + nf( filterRadius, 0, 1 );
    }


    
    fill( 255 );
    noStroke();
    rect( 10-4, winHeight-10-16-1, 16*msg.length()+4+4, 16+6 );

    fill( 0 );
    stroke( 255 );
    text( msg, 10, winHeight-10 );
}
 

// Handle a key press

void keyPressed() {
    switch (key) {

    case CODED:
        switch (keyCode) {
        case RIGHT:
            println( "starting forward FFT" );
            forwardFFT2D();
            showSpectrumImage();
            println( "done forward FFT" );
            break;
        case LEFT:
            println( "starting inverse FFT" );
            inverseFFT2D();
            showSpatialImage();
            println( "done inverse FFT" );
            break;
        }
        break;

    case 'b':
        filterType = (filterType+1) % 3;
        break;
        
    case 'i':
        selectInput( "Select an image", "readNewImage", new File(dataPath("data")) );
        break;

    case '-':
    case '_':
        logFactor /= 1.2;
        showSpectrumImage();
        break;
    case '+':
    case '=':
        logFactor *= 1.2;
        showSpectrumImage();
        break;

    case '>':
    case '.':
        filterRadius *= 1.2;
        break;
    case '<':
    case ',':
        filterRadius /= 1.2;
        break;

    case 'f':
        showSpatialImage();
        break;
    case 'F':
        showSpectrumImage();
        break;

    case '?':
        println( " i   get image from file" );
        println( "-->  compute forward Fourier transform" );
        println( "<--  compute inverse Fourier transform" );
        println( " f   show spatial image (no computation done)" );
        println( " F   show magnitude of Fourier spectrum (no computation done)" );
        println( " -   decrease spectrum enhancement" );
        println( " +   increase spectrum enhancement" );
        println( " b   cycle through blurring filters" );
        println( " <   decrease blurring filter radius" );
        println( " >   increase blurring filter radius" );
        println( " ?   help" );
        break;
    }
}


// Handle a mouse click to modify the frequency signal


void mouseClicked()

{
    if (showingFT) {
        applyBlurFilter();
        showSpectrumImage();
    }
}



// Read an image, resize it to NxN with N a power of two, and store it
// in the spatialSignal array.

void readImage() {
    
    img = loadImage( imageFilename );

    // Find dimensions that are the smallest powers of two at least as
    // large as the image dimensions.

    int nRows = 1;
    while (nRows < img.height)
        nRows = nRows << 1;

    int nCols = 1;
    while (nCols < img.width)
        nCols = nCols << 1;

    // Pick greater of the dimensions.  Store in global var N.

    N = nRows;
    if (N < nCols)
        N = nCols;

    // Create a padded NxN image of the loaded image

    outputImage = createImage( N, N, RGB );

    outputImage.loadPixels();

    for (int i=0; i<N*N; i++)
        outputImage.pixels[i] = color( 0, 0, 0 );

    outputImage.set( (N-img.width)/2, (N-img.height)/2, img );

    outputImage.updatePixels();

    // Keep a copy of the resized image

    img = outputImage.get();
    
    // Copy image intensities into spatialSignal array

    spatialSignal = new float[N*N];

    for (int i=0; i<N*N; i++) 
        spatialSignal[i] = (red( RGB2YCrCb( outputImage.pixels[i] ) ) - minLum) / (float) (maxLum-minLum);

    // Make space for frequency-domain signal (for use later on)

    frequencySignal = new Complex[N*N];

    for (int i=0; i<N*N; i++)
        frequencySignal[i] = new Complex(0,0);
}


void readNewImage( File file ) {
    imageFilename = file.getPath();
    readImage();
    showSpatialImage();
}


// Convert RGB to YCrCb
//
// See the specification "CCIR 601"
//
// Computation:
//       R     G     B         Y       Cr         Cb
//     [0,1] [0,1] [0,1] <-> [0,1] [-0.5,0.5] [-0.5,0.5]
//
// Storage:
//        R       G       B           Y        Cr       Cb
//     [0,255] [0,255] [0,255] <-> [16,235] [16,239] [16,239] 


color RGB2YCrCb( color rgb ) {

    // Get RGB components in [0,1]x[0,1]x[0,1]

    float r = red(   rgb )/255;
    float g = green( rgb )/255;
    float b = blue(  rgb )/255;

    // Compute YCrCb components in [0,1]x[-0.5,0.5]x[0.5,0.5]

    float y  = 0.299*r + 0.587*g + 0.114*b;
    float cr = (b-y)*0.565;
    float cb = (r-y)*0.713;

    // Normalize YUV components into [16,235]x[16,239]x[16,239]
    // (as specified by CCIR 601)
    
    y = y*219+16;
    cr = (cr+0.5)*223+16;
    cb = (cb+0.5)*223+16;

    return color( round(y), round(cr), round(cb) );
}


// Convert YUV to RGB
//
// See the specification "CCIR 601"

color YCrCb2RGB( color ycrcb ) {

    // Get YCrCb components in [0,255]x[0,255]x[0,255]
    
    float y  = red(   ycrcb );
    float cr = green( ycrcb );
    float cb = blue(  ycrcb );

    // Denormalize YCrCb components into [0,1]x[-0.5,0.5]x[0.5,0.5]

    y = (y-16)/219;
    cr = (cr-16)/223-0.5;
    cb = (cb-16)/223-0.5;

    // Compute RGB components in [0,1]x[0,1]x[0,1]

    float r = y               + 1.402520*cb;
    float g = y - 0.343731*cr - 0.714403*cb;
    float b = y + 1.769910*cr;

    // Return RGB in [0,255]x[0,255]x[0,255]

    return color( round(r*255), round(g*255), round(b*255) );
}


// Complex number class

class Complex {

    float real, imag;

    Complex( float r, float i ) {
        real = r;
        imag = i;
    }

    Complex add( Complex c ) {
        return new Complex( real + c.real, imag + c.imag );
    }

    Complex subtract( Complex c ) {
        return new Complex( real - c.real, imag - c.imag );
    }

    Complex mult( Complex c ) {
        return new Complex( real * c.real - imag * c.imag, real * c.imag + imag * c.real );
    }
}


// Bit reversal
// 
// From unifft.c @ MIT

int bitrev( int inp, int numbits )
    
{
    int rev = 0;
    
    for (int i=0; i < numbits; i++) {
        rev = (rev << 1) | (inp & 1);
        inp >>= 1;
    }
    
    return rev;
}


// 1D FFT of an N-sample signal stored in x
//
// N must be a power of two.


float PI = 3.1415926;


void FFT( Complex x[], int N )

{
    // bit-reverse the elements of x
    //
    // From unifft.c @ MIT

    { 
        int log2n = 0;
        int nn = N;
        while (nn > 1) {
            log2n++;
            nn = nn >> 1;
        }

        for (int i=0; i<N; i++) {
            int bri = bitrev(i, log2n);
            if (bri > i) {
                Complex temp = x[i];
                x[i] = x[bri];
                x[bri] = temp;
            }
        }  
    } 

    // YOUR CODE HERE
}


// 2D FFT of an NxN signal stored in x.
//
// x[] is a 1D array of N*N elements, which should be interpreted as a
// 2D array of N rows and N columns.  The element at row r, column c
// is x[c+r*N].
//
// N is a power of two.


void FFT2D( Complex x[], int N )

{
    Complex[] signal = new Complex[N];

    // Compute FFTs of the rows
    // 
    // Copy each row of x into a 1D signal, perform the 1D FFT, then
    // copy the signal back to that row of x.

    // YOUR CODE HERE

    // Compute FFTs of the columns
    //
    // Copy each column of x into a 1D signal, perform the 1D FFT, then
    // copy the signal back to that column of x.

    // YOUR CODE HERE
}



// Compute NxN frequencySignal from NxN spatialSignal.
//
// N, frequencySignal, and spatialSignal are global variables.

void forwardFFT2D()

{
    // Copy spatialSignal into frequencySignal.  Both spatialSignal
    // and frequencySignal are 1D vectors of length NxN representing
    // an NxN image.  Apply transform to centre F.T. at (N/2,N/2).

    // YOUR CODE HERE
 
    // Apply 2D F.T. to frequencySignal

    // YOUR CODE HERE
}


// Compute NxN spatialSignal from NxN frequencySignal.
//
// N, frequencySignal, and spatialSignal are global variables.
//
// Use the complex conjugate to leverage the existing FFT2D routine.

void inverseFFT2D()

{
    // Form complex conjugate of frequencySignal

    // YOUR CODE HERE

    // Apply 2D F.T.

    // YOUR CODE HERE

    // Copy into spatialSignal, normalizing and de-centring

    // YOUR CODE HERE
}


// Draw the FFT spectrum in outputImage

void showSpectrumImage()

{
    // Gather F.T. magnitudes in an array and find max (excluding DC component)

    float[] magnitudes = new float[N*N];

    int i = 0;
    while (i < N*N/2) {
        Complex c = frequencySignal[i];
        float val = sqrt( c.real * c.real +
                              c.imag * c.imag );
        magnitudes[i] = val;
        i++;
    }

    magnitudes[i++] = 0;        // DC component

    while (i < N*N) {
        magnitudes[i] = sqrt( frequencySignal[i].real * frequencySignal[i].real +
                              frequencySignal[i].imag * frequencySignal[i].imag );
        i++;
    }

    // Find max

    float max = 0;

    for (i=0; i<N*N; i++)
        if (magnitudes[i] > max)
            max = magnitudes[i];

    magnitudes[N*N/2] = max;

    // Render spectrum to outputImage

    outputImage.loadPixels();

    for (i=0; i<N*N; i++) {
        int mag = round( 255 * log( logFactor * magnitudes[i]/max + 1 ) / log( logFactor + 1 ));
        outputImage.pixels[i] = color( mag, mag, mag );
    }

    outputImage.updatePixels();

    showingFT = true;
}



// Draw the spatial image
//
// Have to modulate the original intensities with those from the F.T.

void showSpatialImage()

{
    // Render original image, with intensities in spatialSignal

    outputImage.loadPixels();
    img.loadPixels();

    for (int i=0; i<N*N; i++) {
        color ycrcb = RGB2YCrCb( img.pixels[i] );
        outputImage.pixels[i] = YCrCb2RGB( color( spatialSignal[i] * (maxLum-minLum) + minLum, green(ycrcb), blue(ycrcb) ) );
    }

    outputImage.updatePixels();

    showingFT = false;
}



// Apply a radially symmetric blurring filter at mouse position in F.T.
//
// These globals are used:
//
//   filterRadius
//   filterType:  0 = Gaussian, 1 = Butterworth, 2 = box

void applyBlurFilter()

{
    for (int r=-round(filterRadius); r<round(filterRadius); r++)
        for (int c=-round(filterRadius); c<round(filterRadius); c++) {

            // Find offset from F.T. center, which is at (N/2,N/2)

            int x = c+mouseX-N/2;
            int y = r+mouseY-N/2;
            
            if (x>-N/2 && x<N/2 && y>-N/2 && y<N/2) { // inside array?
                
                float dist = sqrt(r*r+c*c)/filterRadius;
                
                if (dist <= 1) { // inside filter radius?

                    // Do four quadrants, symmetrically about (N/2,N/2)

                    for (int i=-1; i<2; i+=2)
                        for (int j=-1; j<2; j+=2) { // apply four times: at +/- x and +/- y.

                            // Apply filter (i.e. multiply filter * F.T.)

                            Complex F = frequencySignal[ (N/2 + i*x) + N*(N/2 + j*y) ];

                            float mag = sqrt( F.real*F.real + F.imag*F.imag );
                            float phase = atan2( F.imag, F.real );

                            // Find filter weight at this distance

                            float weight;

                            switch (filterType) {
                            case 0:
                                weight = exp(-4*dist*dist); // Gaussian
                                break;
                            case 1:
                                weight = 1.0 / pow(1.0 + dist, 2*2); // Butterworth
                                break;
                            case 2:
                                weight = 1; // box
                                break;
                            default:
                                weight = 0;
                                break;
                            }

                            // Update magnitude

                            mag *= (1 - weight);

                            frequencySignal[(N/2 + i*x) + N*(N/2 + j*y)] = new Complex( mag * cos(phase), mag * sin(phase) );
                        }
                }
            }
        }
}
