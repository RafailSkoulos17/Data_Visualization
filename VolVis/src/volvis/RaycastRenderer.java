// Version for students

package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import gui.VolVisApplication;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;
import java.util.*;
import util.TrackballInteractor;
import java.awt.event.*;
import volvis.Renderer;

import java.awt.Color;
import java.awt.*;


/**
 *
 * @author michel
 *  Edit by AVilanova & Nicola Pezzotti
 * 
 * 
 * Main functions to implement the volume rendering
 */


//////////////////////////////////////////////////////////////////////
///////////////// CONTAINS FUNCTIONS TO BE IMPLEMENTED ///////////////
//////////////////////////////////////////////////////////////////////

public class RaycastRenderer extends Renderer implements TFChangeListener {


// attributes
    
    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunction2D tFunc2D;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    Renderer render;
    private boolean mipMode = false;
    private boolean slicerMode = true;
    private boolean compositingMode = false;
    private boolean tf2dMode = false;
    private boolean shadingMode = false;
    private boolean isoMode = false;
    private float iso_value=95; 
    // This is a work around
    private float res_factor = 1.0f;
    private float max_res_factor=0.25f;
    private TFColor isoColor; 
    
    private int prev_mouse_pos = 0;
    private Point current_mouse_pos = new Point();

    
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE MODIFIED    /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    //Function that updates the "image" attribute (result of renderings)
    // using the slicing technique. 
    
    void slicer(double[] viewMatrix) {
        	
        // we start by clearing the image
        resetImage();

        // vector uVec and vVec define the view plane, 
        // perpendicular to the view vector viewVec which is going from the view point towards the object
        // uVec contains the up vector of the camera in world coordinates (image vertical)
        // vVec contains the horizontal vector in world coordinates (image horizontal)
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        getViewPlaneVectors(viewMatrix,viewVec,uVec,vVec);

        // The result of the visualization is saved in an image(texture)
        // we update the vector according to the resolution factor
        // If the resolution is 0.25 we will sample 4 times more points. 
        
        for(int k=0;k<3;k++)
        {
            uVec[k]=res_factor*uVec[k];
            vVec[k]=res_factor*vVec[k];
        }
        
        // compute the volume center
        double[] volumeCenter = new double[3];
        computeVolumeCenter(volumeCenter);

        // Here will be stored the 3D coordinates of every pixel in the plane 
        double[] pixelCoord = new double[3];

        // We get the size of the image/texture we will be puting the result of the 
        // volume rendering operation.
        int imageW=image.getWidth();
        int imageH=image.getHeight();
        

        // Center of the image/texture 
        int[] imageCenter = new int[2];
        imageCenter[0]= imageW/2;
        imageCenter[1]= imageH/2;
        
        // imageW/ image H contains the real width of the image we will use given the resolution. 
        //The resolution is generated once based on the maximum resolution.
        imageW = (int) (imageW*((max_res_factor/res_factor)));
        imageH = (int) (imageH*((max_res_factor/res_factor)));

        
        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();

        // Color that will be used as a result 
        TFColor pixelColor = new TFColor();
        
        // Auxiliar color
        TFColor colorAux;

        // Contains the voxel value of interest
        int val;
        
        //Iterate on every pixel
        for (int j = imageCenter[1] - imageH/2; j < imageCenter[1] + imageH/2; j++) {
            for (int i =  imageCenter[0] - imageW/2; i <imageCenter[0] + imageW/2; i++) {
                
                // computes the pixelCoord which contains the 3D coordinates of the pixels (i,j)
                computePixelCoordinatesFloat(pixelCoord,volumeCenter,uVec,vVec,i,j);
                
                //we now have to get the value for the in the 3D volume for the pixel
                //we can use a nearest neighbor implementation like this:
//                val = volume.getVoxelNN(pixelCoord);
                
                //you have also the function getVoxelLinearInterpolated in Volume.java          
//                val = (int) volume.getVoxelLinearInterpolate(pixelCoord);
                    
                
                //you have to implement this function below to get the cubic interpolation
                val = (int) volume.getVoxelTriCubicInterpolate(pixelCoord);
                
                
                // Map the intensity to a grey value by linear scaling
                 pixelColor.r = (val/max);
                 pixelColor.g = pixelColor.r;
                 pixelColor.b = pixelColor.r;

                // the following instruction makes intensity 0 completely transparent and the rest opaque
                // pixelColor.a = val > 0 ? 1.0 : 0.0;   
                
                // Alternatively, apply the transfer function to obtain a color using the tFunc attribute
//                 colorAux= tFunc.getColor(val);
//                 
//                 pixelColor.r=colorAux.r;
//                 pixelColor.g=colorAux.g;
//                 pixelColor.b=colorAux.b;
//                 pixelColor.a=colorAux.a; 
                 
                 
                // IMPORTANT: You can also simply use pixelColor = tFunc.getColor(val); However then you copy by reference and this means that if you change 
                // pixelColor you will be actually changing the transfer function So BE CAREFUL when you do this kind of assignments

                //BufferedImage/image/texture expects a pixel color packed as ARGB in an int
                //use the function computeImageColor to convert your double color in the range 0-1 to the format need by the image
                int pixelColor_i = computeImageColor(pixelColor.r,pixelColor.g,pixelColor.b,pixelColor.a);
    
                image.setRGB(i, j, pixelColor_i);
            }
        }
}
    

    
    //Do NOT modify this function
    //
    //Function that updates the "image" attribute using the MIP raycasting
    //It returns the color assigned to a ray/pixel given it's starting point (entryPoint) and the direction of the ray(rayVector).
    // exitPoint is the last point.
    //ray must be sampled with a distance defined by the sampleStep
   
    int traceRayMIP(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {
        
    	//compute the increment and the number of samples
        double[] increments = new double[3];
        VectorMath.setVector(increments, rayVector[0] * sampleStep, rayVector[1] * sampleStep, rayVector[2] * sampleStep);
        
        // Compute the number of times we need to sample
        double distance = VectorMath.distance(entryPoint, exitPoint);
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);

        //the current position is initialized as the entry point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, entryPoint[0], entryPoint[1], entryPoint[2]);
       
        double maximum = 0;
        do {
            double value = volume.getVoxelLinearInterpolate(currentPos)/255.; 
            if (value > maximum) {
                maximum = value;
            }
            for (int i = 0; i < 3; i++) {
                currentPos[i] += increments[i];
            }
            nrSamples--;
        } while (nrSamples > 0);

        double alpha;
        double r, g, b;
        if (maximum > 0.0) { // if the maximum = 0 make the voxel transparent
            alpha = 1.0;
        } else {
            alpha = 0.0;
        }
        r = g = b = maximum;
        int color = computeImageColor(r,g,b,alpha);
        return color;
    }
    
          
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    //Function that updates the "image" attribute using the Isosurface raycasting
    //It returns the color assigned to a ray/pixel given it's starting point (entryPoint) and the direction of the ray(rayVector).
    // exitPoint is the last point.
    //ray must be sampled with a distance defined by the sampleStep
   
   int traceRayIso(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {
       
        TFColor shading = new TFColor();
        TFColor newcolor = new TFColor();
        
        VoxelGradient gradient = new VoxelGradient();
        
        double[] viewVector = new double[3];
        VectorMath.setVector(viewVector, -rayVector[0], -rayVector[1], -rayVector[2]);
        double[] lightVector = new double[3];
        //We define the light vector as directed toward the view point (which is the source of the light)
        // another light vector would be possible
         VectorMath.setVector(lightVector, rayVector[0], rayVector[1], rayVector[2]);
         
        // defines if bisection will be used
        boolean bisection = false;
        //defines if there was a hit
        boolean found=false;
        
        //the current position is initialized as the entry point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, entryPoint[0], entryPoint[1], entryPoint[2]);
        double[] increments = new double[3];
        VectorMath.setVector(increments, rayVector[0] * sampleStep, rayVector[1] * sampleStep, rayVector[2] * sampleStep);
        double[] prevPos = new double[3];
        VectorMath.setVector(prevPos, currentPos[0], currentPos[1], currentPos[2]);
        double[] best_position = new double[3];
        
        float previousvalue = 0;
        float value = 0;
        
       // compute the number of times we need to sample
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);
        
        //sample through the ray
        while (!found && nrSamples>0){
            previousvalue = value;
            //get the current intensity value
            value = volume.getVoxelLinearInterpolate(currentPos); 
            if (value >= iso_value) { //if there was a hit
                found=true; 
                if (bisection){ //if we use bisection mode
                    best_position = bisection_accuracy (currentPos, prevPos, previousvalue, value, iso_value);
                    gradient = gradients.getGradient(best_position);
                }
                else{
                    gradient = gradients.getGradient(currentPos); 
                }
                break;
            }
            // set the previous position
            VectorMath.setVector(prevPos, currentPos[0], currentPos[1], currentPos[2]);
            
            // move to the next position
            for (int i = 0; i < 3; i++) {
                currentPos[i] += increments[i];
            }
            nrSamples--;
        }

        //Initialization of the colors as floating point values
        double alpha;
        
        // if there was a miss, make the voxel transarent
        if (found) { 
            alpha = 1.0;
        } else {
            alpha = 0.0;
        }             
        
        //isoColor contains the isosurface color from the interface
        newcolor.r = isoColor.r;
        newcolor.g = isoColor.g;
        newcolor.b =  isoColor.b;
        
        //if shading mode is enabled
        if(shadingMode){
                //if the voxel is not transparent we apply shading using Phong's model
                if(alpha==1.0){
                    //calculate the color and the opacity after the application of Phong Shading
                    shading=computePhongShading(newcolor, gradient, lightVector, viewVector);
                    
                    newcolor.r=shading.r;
                    newcolor.g=shading.g;
                    newcolor.b=shading.b;
                    newcolor.a=shading.a;
                    alpha = shading.a;

                }
            }
        
        //compute the color
        int color = computeImageColor(newcolor.r,newcolor.g,newcolor.b,alpha);
        return color;
    }
   
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    
   // Given the current sample position, increment vector of the sample (vector from previous sample to current sample) and sample Step. 
   // Previous sample value and current sample value, isovalue value
   // The function should search for a position where the iso_value passes that it is more precise.
   // Some unused argumenets were removed
   double[] bisection_accuracy (double[] currentPos, double[] prevPos, float previousvalue,float value, float iso_value)
   {    
       
        double[] middlePos = new double[3];
        // define the threshold we want for the difference between iso_value and the current value
        double isoTol = 1;
        // define maximum number of bisections
        
        int numBisections = 20;
        int middleValue = 0;
       
       for (int i=0; i<numBisections; i++){
            
           // get the middle point between current and previous position
            VectorMath.setVector(middlePos, (prevPos[0] + currentPos[0])/2.0, (prevPos[1] + currentPos[1])/2.0, (prevPos[2] + currentPos[2])/2.0);
            // get the intensity for the middle point
            middleValue = (int) volume.getVoxelLinearInterpolate(middlePos);
            
            // compares the isovalue and the value in the middle point
            // and it devided the interval in the right way
            if (Math.abs(middleValue - iso_value) <= isoTol){ //if we are close to the isovalue return this position
                return middlePos;
            }
            else if (middleValue > iso_value){ // if the value in the middle point is bigger than isovalue, 
                if (previousvalue > value){ // if the intensity is decreasing
                    prevPos = middlePos;
                }
                else{
                    currentPos = middlePos;
                }
            }
            else if (middleValue < iso_value){ // if the value in the middle point is smaller than isovalue
                if (previousvalue > value){ // if the intensity is decreasing
                    currentPos = middlePos; 
                }
                else{
                     prevPos = middlePos;     
                }
            }     
        }
        return middlePos;
   }
    
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    //Function that updates the "image" attribute using the compositing// accumulatted raycasting
    //It returns the color assigned to a ray/pixel given it's starting point (entryPoint) and the direction of the ray(rayVector).
    // exitPoint is the last point.
    //ray must be sampled with a distance defined by the sampleStep
   
    TFColor BackToFront(TFColor AccumulatedColor, TFColor newColor, double alpha){
        AccumulatedColor.r= newColor.r * alpha + (1.0 - alpha) * AccumulatedColor.r;
        AccumulatedColor.g= newColor.g * alpha + (1.0 - alpha) * AccumulatedColor.g;
        AccumulatedColor.b= newColor.b * alpha + (1.0 - alpha) * AccumulatedColor.b;
        AccumulatedColor.a = alpha + (1.0 - alpha) * AccumulatedColor.a;
     return AccumulatedColor;
   }
   
    int traceRayComposite(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {
        
        TFColor voxel_color = new TFColor();
        TFColor newColor=new TFColor();
        TFColor shading=new TFColor();
        
        //initialize variable color where we store the result of color composition
        TFColor AccumulatedColor=new TFColor();
        
        //initialize the light vector
        double[] lightVector = new double[3];
        VectorMath.setVector(lightVector, rayVector[0], rayVector[1], rayVector[2]);
        
        //initialize the view vector
        double[] viewVector = new double[3];
        VectorMath.setVector(viewVector, -rayVector[0], -rayVector[1], -rayVector[2]);

        //compute the increment and the number of samples
        double[] increments = new double[3];
        VectorMath.setVector(increments, rayVector[0] * sampleStep, rayVector[1] * sampleStep, rayVector[2] * sampleStep);
        
        
        //the current position is initialized as the entry point
        double[] currentPos = new double[3];
        double voxel_intensity;
        double alpha = 0; 
        
        // Define the distance and direction between each sample step along the three axes
    	double stepDistX = -rayVector [0] * sampleStep;
        double stepDistY = -rayVector [1] * sampleStep;
        double stepDistZ = -rayVector [2] * sampleStep;

        // Define array that holds the 3 coordinates of each sample.
        double[] voxelCoord = new double[3];
        //define the first sample which will be the last one (due to back to front color composing)
        voxelCoord[0] = exitPoint[0];
        voxelCoord[1] = exitPoint[1];
        voxelCoord[2] = exitPoint[2];
        
        double[] prevCoord = new double[3];
        double[] nexCoord = new double[3];
        
        VectorMath.setVector(currentPos, exitPoint[0], exitPoint[1], exitPoint[2]);
        int resultColor=0; 

        // Compute the number of times we need to sample
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);
        
        //initialize all colors and opacity with zero
        AccumulatedColor.r=0.0;
        AccumulatedColor.b=0.0;
        AccumulatedColor.g=0.0;
        AccumulatedColor.a=0.0;
        
        

        
        //the light vector is directed toward the view point (which is the source of the light)
        // itarate through the ray
        do {

            // get the gradient and the intensity of the voxel
            voxel_intensity = volume.getVoxelLinearInterpolate(currentPos);
            VoxelGradient gradient = gradients.getGradient(currentPos);

            if (compositingMode){
                //based on the sample intensity, we get its color and opacity through the transfer function
                voxel_color = tFunc.getColor((int)voxel_intensity);
                alpha = voxel_color.a;
                //assign the color to our variable newColor
                newColor.r =voxel_color.r;
                newColor.g =voxel_color.g;
                newColor.b =voxel_color.b;
                newColor.a =voxel_color.a;
            }
            
            //if 2d transfer function mode is enabled
            if (tf2dMode){
                
                //get the color of the voxel via the 2d transfer function
                voxel_color=tFunc2D.color;
                
                 //assign the color to our variable newColor
                newColor.a = voxel_color.a;
                newColor.r = voxel_color.r;
                newColor.g = voxel_color.g;
                newColor.b = voxel_color.b;
                
                //for this mode we have to calculate the opacity from the triangular classification widget
                //the following function gets the value of the material, the factor r, the voxel value and the gradient magnitude
                // and returns the opacity to be assgned
                alpha = computeOpacity2DTF(tFunc2D.baseIntensity,tFunc2D.radius, voxel_intensity, gradient.mag);
                newColor.a = alpha;
            }
            
            //if shading mode is enabled
            if(shadingMode){
                
                //if the voxel is not transparent we apply shading using Phong's shading model
                if(alpha>0.0){
                    //calculate the color and opacity after the application of phong shading
                    shading=computePhongShading(newColor, gradient, lightVector, viewVector);
                    
                    //assign the color to our variable newColor
                    newColor.r=shading.r;
                    newColor.g=shading.g;
                    newColor.b=shading.b;
                    newColor.a=shading.a;
                    alpha = shading.a;

                }
            }
            
            // apply back to front compositing
            AccumulatedColor=BackToFront(AccumulatedColor,newColor,alpha);
           
            //implementation of early ray termination to speed up our calculations
            if (newColor.a>=0.98 && newColor.a<=1.0){
            	//if opacity of sample is approximately 1.0 we stop calculating samples along the ray
                break;
            }
            //update current position
            for (int i = 0; i < 3; i++) {
                currentPos[i] -= increments[i]; 
            }
                        
            nrSamples--;
        } while (nrSamples > 0);
        
        AccumulatedColor.a = 1.0;
       //compute the final color to be assigned
       resultColor=computeImageColor(AccumulatedColor.r,AccumulatedColor.g,AccumulatedColor.b,AccumulatedColor.a);
       return resultColor;
    }
    

    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // Compute Phong Shading given the voxel color (material color), the gradient, the light vector and view vector 
    // We assume that there is a typo in the arguments of the function, that is, instead of rayVector should be viewVector
    TFColor computePhongShading(TFColor voxel_color, VoxelGradient gradient, double[] lightVector,
         double[] viewVector) {
        
        // initialize the shading color for the given voxel
	TFColor shadingColor=new TFColor(voxel_color.r,voxel_color.g,voxel_color.b,voxel_color.a);
        
        //parameters as given by the exercise
        double kd=0.7;
	double ka=0.1;
	double ks=0.2;
	double n_power=100;

        //normal vector 
	double[] normal= new double[3];
        //normalize the gradient by dividing it by its magnitude and keep it to the normal vector
	VectorMath.setVector(normal,gradient.x/gradient.mag, gradient.y/gradient.mag, gradient.z/gradient.mag);

	//calculate ambient and diffuse
        
        //calculate the cosine of the angle formed by lightVector and normal vector
	double diffuse=VectorMath.dotproduct(normal,lightVector);

        // calculate the color that results from ambient and diffuse based on the equation Ld*kd*ci,d(L.n)+La*ka*ci,a
        // where we assumed that La,Ld are equal to one
	if (diffuse > 0){
            shadingColor.r = shadingColor.r * diffuse * kd + voxel_color.r * ka ;
            shadingColor.g = shadingColor.g * diffuse * kd + voxel_color.g * ka;
            shadingColor.b = shadingColor.b * diffuse * kd + voxel_color.b * ka;
	}

	//calculate the r vector based on the equation  R=2(n.L)*n-L
      	double [] rVector=new double[3];
        for(int i=0;i<3;i++){
            rVector[i]=2*diffuse*normal[i]-lightVector[i];
        }

	//Calculate specular and VR product(this product is equal to the cosine of the angle between the viewVector and r vector)
	double specular=VectorMath.dotproduct(normal,rVector);
        double vrProduct=VectorMath.dotproduct(viewVector,rVector);

        //Calculate specular based on the equation Ls*ks*ci,s*(V.R)^n
	if (specular > 0) {
            //following is the procedure for specular based on Phong's model
//            shadingColor.r += ks * Math.pow(vrProduct, n_power);
//            shadingColor.g += ks * Math.pow(vrProduct, n_power);
//            shadingColor.b += ks * Math.pow(vrProduct, n_power);
            
            //we tried implemented an alternative for Phong's Model which is faster
            shadingColor.r += ks * ((vrProduct/(n_power-n_power*vrProduct +vrProduct)));
            shadingColor.g += ks * ((vrProduct/(n_power-n_power*vrProduct +vrProduct)));
            shadingColor.b += ks * ((vrProduct/(n_power-n_power*vrProduct +vrProduct)));
            
            // we used the absolute values, because of some negative values which were observed
            shadingColor.r = Math.abs(shadingColor.r);
            shadingColor.g = Math.abs(shadingColor.g);
            shadingColor.b = Math.abs(shadingColor.b);

        }

        shadingColor.r = shadingColor.r > 1.0 ? 1.0 : shadingColor.r;
        shadingColor.g = shadingColor.g > 1.0 ? 1.0 : shadingColor.g;
        shadingColor.b = shadingColor.b > 1.0 ? 1.0 : shadingColor.b;
        
        return shadingColor;
    }
    
    
    //////////////////////////////////////////////////////////////////////
    ///////////////// LIMITED MODIFICATION IS NEEDED /////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // Implements the basic tracing of rays trough the image and given the
    // camera transformation
    // It calls the functions depending on the raycasting mode
  
    void raycast(double[] viewMatrix) {
        
    	//data allocation
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] pixelCoord = new double[3];
        double[] entryPoint = new double[3];
        double[] exitPoint = new double[3];
        
        // increment in the pixel domain in pixel units
        int increment = 1;
        // sample step in voxel units
        int sampleStep = 1;
        
         //We increases the sample step to 5 so we get less samples per ray 
        //and the number of pixels we look at a time to 2, when user is interacting with the simulator
        //This way the simulator is more responsive and does not lag
        if (this.interactiveMode){
            increment = 2;
            sampleStep = 5;
        }
        
        // reset the image to black
        resetImage();
        

        // vector uVec and vVec define the view plane, 
        // perpendicular to the view vector viewVec which is going from the view point towards the object
        // uVec contains the up vector of the camera in world coordinates (image vertical)
        // vVec contains the horizontal vector in world coordinates (image horizontal)
        getViewPlaneVectors(viewMatrix,viewVec,uVec,vVec);
        
        
        // The result of the visualization is saved in an image(texture)
        // we update the vector according to the resolution factor
        // If the resolution is 0.25 we will sample 4 times more points. 
        for(int k=0;k<3;k++)
        {
            uVec[k]=res_factor*uVec[k];
            vVec[k]=res_factor*vVec[k];
        }
        
        
        
       // We get the size of the image/texture we will be puting the result of the 
        // volume rendering operation.
        int imageW=image.getWidth();
        int imageH=image.getHeight();

        int[] imageCenter = new int[2];
        // Center of the image/texture 
        imageCenter[0]= imageW/2;
        imageCenter[1]= imageH/2;
        
        // imageW/ image H contains the real width of the image we will use given the resolution. 
        //The resolution is generated once based on the maximum resolution.
        imageW = (int) (imageW*((max_res_factor/res_factor)));
        imageH = (int) (imageH*((max_res_factor/res_factor)));
    
        //The rayVector is pointing towards the scene
        double[] rayVector = new double[3];
        rayVector[0]=-viewVec[0];
        rayVector[1]=-viewVec[1];
        rayVector[2]=-viewVec[2];

             
        // compute the volume center
        double[] volumeCenter = new double[3];
        computeVolumeCenter(volumeCenter);

        
        // ray computation for each pixel
        for (int j = imageCenter[1] - imageH/2; j < imageCenter[1] + imageH/2; j += increment) {
            for (int i =  imageCenter[0] - imageW/2; i <imageCenter[0] + imageW/2; i += increment) {
                // compute starting points of rays in a plane shifted backwards to a position behind the data set
            	computePixelCoordinatesBehindFloat(pixelCoord,viewVec,uVec,vVec,i,j);
            	// compute the entry and exit point of the ray
                computeEntryAndExit(pixelCoord, rayVector, entryPoint, exitPoint);
                if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                    int val = 0;
                    if (compositingMode || tf2dMode) {
                        val = traceRayComposite(entryPoint, exitPoint, rayVector, sampleStep);
                    } else if (mipMode) {
                        val = traceRayMIP(entryPoint, exitPoint, rayVector, sampleStep);
                    } else if (isoMode){
                        val= traceRayIso(entryPoint,exitPoint,rayVector, sampleStep);
                    }
                    for (int ii = i; ii < i + increment; ii++) {
                        for (int jj = j; jj < j + increment; jj++) {
                            image.setRGB(ii, jj, val);
                        }
                    }
                }

            }
        }
    }


//////////////////////////////////////////////////////////////////////
///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
////////////////////////////////////////////////////////////////////// 
// Compute the opacity based on the value of the pixel and the values of the
// triangle widget tFunc2D contains the values of the baseintensity and radius
// tFunc2D.baseIntensity, tFunc2D.radius they are in image intensity units

 
    
    
public double computeOpacity2DTF(double material_value, double material_r,
        double voxelValue, double gradMagnitude) {
    
    //initialize opacity to zero
    double opacity = 0.0;
    
    //get the coordinbates of the triangle's vertexes a, b and c
    double[] a = {material_value - material_r, 173.1};
    double[] b = {material_value + material_r, 173.1};
    double[] c = {material_value, 0};
    
    //get the coordinates of the center of mass(centroid) of the triangle
    double x_center = (a[0] + b[0] + c[0]) / 3;
    double y_center = (a[1] + b[1] + c[1]) / 3;
    double[] center = {x_center, y_center};

    
    //get the coordinates of the current point in the grid
    double[] currPoint = {voxelValue, gradMagnitude};

    // if the current point is not inside the triangle
    if (!isInside(a[0],a[1], b[0], b[1], c[0], c[1], currPoint[0], currPoint[1])){
        opacity = 0;
    }
    else{ // if the current point is  inside the triangle
          double y = currPoint[1];
          double x = currPoint[0];
          if (x == material_value){ // if voxel value is equal to the material value
            opacity = 1;
          }
          else if (x > material_value){ // if it is on the right part of the triangle
              
              // calculate the distance form the marerial value and the maximum value of this distance
              double max_dist =  calculateDistanceFromLineXaxis(b[0], b[1], c[0], c[1], y, material_value);
              double distance  = Math.abs(x - material_value);
        
    //        opacity is linearly decreasing as the distance from the material value increases
              opacity = 1 - (distance/max_dist);
          }
          else{ // if it is on the left part of the triangle
              
              // calculate the distance form the marerial value and the maximum value of this distance
              double max_dist =  calculateDistanceFromLineXaxis(a[0], a[1], c[0], c[1], y, material_value);
              double distance  = Math.abs(x - material_value);
              
    //        opacity is linearly decreasing as the distance from the material value increases
              opacity = 1 - (distance/max_dist);
          }
        }
    
    return opacity;
    }  


    // a utility function to calculate area of triangle  
    //formed by (x1, y1) (x2, y2) and (x3, y3)
    static double area(double x1, double y1, double x2, double y2, 
                                        double x3, double y3) 
    { 
       return Math.abs((x1*(y2-y3) + x2*(y3-y1)+ 
                                    x3*(y1-y2))/2.0); 
    } 
       
    // a function to check whether point P(x, y) lies 
    // inside the triangle formed by A(x1, y1), 
    // B(x2, y2) and C(x3, y3) */
    static boolean isInside(double x1, double y1, double x2, 
                double y2, double x3, double y3, double x, double y) 
    {    
       // calculate area of triangle ABC
        double A = area (x1, y1, x2, y2, x3, y3); 
       
        // calculate area of triangle PBC  
        double A1 = area (x, y, x2, y2, x3, y3); 
       
        // calculate area of triangle PAC
        double A2 = area (x1, y1, x, y, x3, y3); 
       
        // calculate area of triangle PAB */   
        double A3 = area (x1, y1, x2, y2, x, y); 
         
        // check if sum of A1, A2 and A3 is same as A
        return (A == A1 + A2 + A3); 
    }
    
    // a function to calculate the Euclidean distance between two points
    public static double calculateDistance(double[] A, double[] B)
    {
        double sum = 0;
        for (int i=0;i<A.length;i++){
            sum += Math.pow((A[i] - B[i]), 2.0); 
        }
    return Math.sqrt(sum);
    }
    
    // a function that calculates the distance  between vixel value and the 
    // maximum value the intensity can have inside the triangle for the given 
    // gradient of the voxel
    public static double calculateDistanceFromLineXaxis(double x1, double y1, double x2,double y2, double y, double material_value){
        double x = (y*(x2-x1) + x1*y2 - y1*x2)/ (y2-y1);
        return Math.abs(x - material_value);
    
    }
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  
    
    //Do NOT modify this function
    int computeImageColor(double r, double g, double b, double a){
		int c_alpha = 	a <= 1.0 ? (int) Math.floor(a * 255) : 255;
        int c_red = 	r <= 1.0 ? (int) Math.floor(r * 255) : 255;
        int c_green = 	g <= 1.0 ? (int) Math.floor(g * 255) : 255;
        int c_blue = 	b <= 1.0 ? (int) Math.floor(b * 255) : 255;
        int pixelColor = getColorInteger(c_red,c_green,c_blue,c_alpha);
        return pixelColor;
	}
    //Do NOT modify this function    
    public void resetImage(){
    	for (int j = 0; j < image.getHeight(); j++) {
	        for (int i = 0; i < image.getWidth(); i++) {
	            image.setRGB(i, j, 0);
	        }
	    }
    }
   
    //Do NOT modify this function
    void computeIncrementsB2F(double[] increments, double[] rayVector, double sampleStep) {
        // we compute a back to front compositing so we start increments in the oposite direction than the pixel ray
    	VectorMath.setVector(increments, -rayVector[0] * sampleStep, -rayVector[1] * sampleStep, -rayVector[2] * sampleStep);
    }
    
    //used by the slicer
    //Do NOT modify this function
    void getViewPlaneVectors(double[] viewMatrix, double viewVec[], double uVec[], double vVec[]) {
            VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
	    VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
	    VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
	}
    
    //used by the slicer	
    //Do NOT modify this function
    void computeVolumeCenter(double volumeCenter[]) {
	VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
    }
	
    //used by the slicer
    //Do NOT modify this function
    void computePixelCoordinatesFloat(double pixelCoord[], double volumeCenter[], double uVec[], double vVec[], float i, float j) {
        // Coordinates of a plane centered at the center of the volume (volumeCenter and oriented according to the plane defined by uVec and vVec
            float imageCenter = image.getWidth()/2;
            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];
    }
        
    //Do NOT modify this function
    void computePixelCoordinates(double pixelCoord[], double volumeCenter[], double uVec[], double vVec[], int i, int j) {
        // Coordinates of a plane centered at the center of the volume (volumeCenter and oriented according to the plane defined by uVec and vVec
            int imageCenter = image.getWidth()/2;
            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];
    }
    //Do NOT modify this function
    void computePixelCoordinatesBehindFloat(double pixelCoord[], double viewVec[], double uVec[], double vVec[], float i, float j) {
            int imageCenter = image.getWidth()/2;
            // Pixel coordinate is calculate having the center (0,0) of the view plane aligned with the center of the volume and moved a distance equivalent
            // to the diaganal to make sure I am far away enough.

            double diagonal = Math.sqrt((volume.getDimX()*volume.getDimX())+(volume.getDimY()*volume.getDimY())+ (volume.getDimZ()*volume.getDimZ()))/2;               
            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * diagonal + volume.getDimX() / 2.0;
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * diagonal + volume.getDimY() / 2.0;
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * diagonal + volume.getDimZ() / 2.0;
    }
    //Do NOT modify this function
    void computePixelCoordinatesBehind(double pixelCoord[], double viewVec[], double uVec[], double vVec[], int i, int j) {
            int imageCenter = image.getWidth()/2;
            // Pixel coordinate is calculate having the center (0,0) of the view plane aligned with the center of the volume and moved a distance equivalent
            // to the diaganal to make sure I am far away enough.

            double diagonal = Math.sqrt((volume.getDimX()*volume.getDimX())+(volume.getDimY()*volume.getDimY())+ (volume.getDimZ()*volume.getDimZ()))/2;               
            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * diagonal + volume.getDimX() / 2.0;
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * diagonal + volume.getDimY() / 2.0;
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * diagonal + volume.getDimZ() / 2.0;
    }
	
    //Do NOT modify this function
    public int getColorInteger(int c_red, int c_green, int c_blue, int c_alpha) {
    	int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
    	return pixelColor;
    } 
    //Do NOT modify this function
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
        isoColor = new TFColor();
        isoColor.r=1.0;isoColor.g=1.0;isoColor.b=0.0;isoColor.a=1.0;
    }
    
     
    //Do NOT modify this function
    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ())* (1/max_res_factor));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
      
                  
           
        // Initialize transferfunction 
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        tFunc.setTestFunc();
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tFunc2D= new TransferFunction2D((short) (volume.getMaximum() / 2), 0.2*volume.getMaximum());
        tfEditor2D = new TransferFunction2DEditor(tFunc2D,volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }
    //Do NOT modify this function
    public RaycastRendererPanel getPanel() {
        return panel;
    }

    //Do NOT modify this function
    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    //Do NOT modify this function
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
    //Do NOT modify this function
    public void setShadingMode(boolean mode) {
        shadingMode = mode;
        changed();
    }
    //Do NOT modify this function
    public void setMIPMode() {
        setMode(false, true, false, false,false);
    }
    //Do NOT modify this function
    public void setSlicerMode() {
        setMode(true, false, false, false,false);
    }
    //Do NOT modify this function
    public void setCompositingMode() {
        setMode(false, false, true, false,false);
    }
    //Do NOT modify this function
    public void setTF2DMode() {
        setMode(false, false, false, true, false);
    }
    //Do NOT modify this function
    public void setIsoSurfaceMode(){
        setMode(false, false, false, false, true);
     }
    //Do NOT modify this function
    public void setIsoValue(float pIsoValue){
         iso_value = pIsoValue;
         if (isoMode){
             changed();
         }
             
     }
    //Do NOT modify this function
    public void setResFactor(int value) {
         float newRes= 1.0f/value;
         if (res_factor != newRes)
         {
             res_factor=newRes;
             if (volume != null) changed();
         }
     }
     
   //Do NOT modify this function
   public void setIsoColor(TFColor newColor)
     {
         this.isoColor.r=newColor.r;
         this.isoColor.g=newColor.g;
         this.isoColor.b=newColor.b;
         if ((volume!=null) && (this.isoMode)) changed();
     }
     
    //Do NOT modify this function
     public float getIsoValue(){
         return iso_value;
     }
    //Do NOT modify this function
    private void setMode(boolean slicer, boolean mip, boolean composite, boolean tf2d, boolean iso) {
        slicerMode = slicer;
        mipMode = mip;
        compositingMode = composite;
        tf2dMode = tf2d;        
        isoMode = iso;
        changed();
    }
    //Do NOT modify this function
    private boolean intersectLinePlane(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection) {

        double[] tmp = new double[3];

        for (int i = 0; i < 3; i++) {
            tmp[i] = plane_pos[i] - line_pos[i];
        }

        double denom = VectorMath.dotproduct(line_dir, plane_normal);
        if (Math.abs(denom) < 1.0e-8) {
            return false;
        }

        double t = VectorMath.dotproduct(tmp, plane_normal) / denom;

        for (int i = 0; i < 3; i++) {
            intersection[i] = line_pos[i] + t * line_dir[i];
        }

        return true;
    }
    //Do NOT modify this function
    private boolean validIntersection(double[] intersection, double xb, double xe, double yb,
            double ye, double zb, double ze) {

        return (((xb - 0.5) <= intersection[0]) && (intersection[0] <= (xe + 0.5))
                && ((yb - 0.5) <= intersection[1]) && (intersection[1] <= (ye + 0.5))
                && ((zb - 0.5) <= intersection[2]) && (intersection[2] <= (ze + 0.5)));

    }
    //Do NOT modify this function
    private void intersectFace(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection,
            double[] entryPoint, double[] exitPoint) {

        boolean intersect = intersectLinePlane(plane_pos, plane_normal, line_pos, line_dir,
                intersection);
        if (intersect) {

            double xpos0 = 0;
            double xpos1 = volume.getDimX();
            double ypos0 = 0;
            double ypos1 = volume.getDimY();
            double zpos0 = 0;
            double zpos1 = volume.getDimZ();

            if (validIntersection(intersection, xpos0, xpos1, ypos0, ypos1,
                    zpos0, zpos1)) {
                if (VectorMath.dotproduct(line_dir, plane_normal) < 0) {
                    entryPoint[0] = intersection[0];
                    entryPoint[1] = intersection[1];
                    entryPoint[2] = intersection[2];
                } else {
                    exitPoint[0] = intersection[0];
                    exitPoint[1] = intersection[1];
                    exitPoint[2] = intersection[2];
                }
            }
        }
    }
    
     
    
    //Do NOT modify this function
    void computeEntryAndExit(double[] p, double[] viewVec, double[] entryPoint, double[] exitPoint) {

        for (int i = 0; i < 3; i++) {
            entryPoint[i] = -1;
            exitPoint[i] = -1;
        }

        double[] plane_pos = new double[3];
        double[] plane_normal = new double[3];
        double[] intersection = new double[3];

        VectorMath.setVector(plane_pos, volume.getDimX(), 0, 0);
        VectorMath.setVector(plane_normal, 1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, -1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, volume.getDimY(), 0);
        VectorMath.setVector(plane_normal, 0, 1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, -1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, volume.getDimZ());
        VectorMath.setVector(plane_normal, 0, 0, 1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, 0, -1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

    }

    //Do NOT modify this function
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();


    }
    //Do NOT modify this function
    @Override
    public void visualize(GL2 gl) {

        double[] viewMatrix = new double[4 * 4];
        
        if (volume == null) {
            return;
        }
        	
         drawBoundingBox(gl);

        
     //    gl.glGetDoublev(GL2.GL_PROJECTION_MATRIX,viewMatrix,0);

         gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);
                      
         
        long startTime = System.currentTimeMillis();
        if (slicerMode) {
            slicer(viewMatrix);    
        } else {
            raycast(viewMatrix);
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);
        
        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = res_factor*image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST);
        gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST);
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(texture.getImageTexCoords().left(), texture.getImageTexCoords().top());
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(texture.getImageTexCoords().left(), texture.getImageTexCoords().bottom());
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(texture.getImageTexCoords().right(), texture.getImageTexCoords().bottom());
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(texture.getImageTexCoords().right(), texture.getImageTexCoords().top());
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }


    }
    
    //Do NOT modify this function
    private BufferedImage image;

    //Do NOT modify this function
    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
    

}