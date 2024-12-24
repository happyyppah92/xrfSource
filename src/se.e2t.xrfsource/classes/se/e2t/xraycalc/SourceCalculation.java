/*
 * SourceCalculation.java
 *
 * Copyright 2019 e2t AB
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the Software
 * is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package se.e2t.xraycalc;

import se.e2t.abscoeffcalculate.AbsCoefficient;
import se.e2t.abscoeffcalculate.Mucal;
import java.util.List;
import se.e2t.xrfsource.spectrumclasses.SpectrumPart;
import se.e2t.xrfsource.spectrumclasses.XraySpectrum;

/**
 *
 * @author Kent Ericsson, e2t AB
 * 
 * This is the super class for calculation of tube spectrum intensities.
 * This class is extended by the classes that actually calculate the intensities.
 */
public abstract class SourceCalculation {

    protected static final double MINIMUM_SLICE = 0.001d;

    public SourceCalculation() {
    }

    public XraySpectrum calculate(Inparameters inParameters) {

        XraySpectrum outputData = new XraySpectrum();

        // Find center wavelength of first full spectrum slice     
        //add by cx
        double endTubeVoltage= inParameters.getTubeVoltage();
        double deltaTV= inParameters.getContinuumIntervalSize();
        double minTubeVoltage=inParameters.getMintubevoltage();
        double centerTV=minTubeVoltage+deltaTV/2.0d;
        double TVuppper=minTubeVoltage+deltaTV;

        //改成输入的是kv不是波长
        if(minTubeVoltage < TVuppper)
        {
            do {
                double intensity= getContiniumIntensity(inParameters, centerTV, deltaTV);
                outputData.addContiniumSlice(
                        new SpectrumPart(centerTV, deltaTV,intensity ));
                centerTV += deltaTV;
                TVuppper += deltaTV;
            } while (TVuppper < endTubeVoltage + deltaTV); // to pass max
        
        }
        // Calculate tube line intensities
        calculateTubeLineIntensities(inParameters, outputData);

        // Adjust intensities depending on window and filter attenuation
        windowFilterAdjustment(inParameters, outputData);

        // Normalize calculated intensities
        //normalizeIntensities(outputData, 1.0d);

        return outputData;
    }
    
    // This method is implemented by the classes extending this class
    /**
     * This method which is implemented by the classes extending this class
     * calculates continuum intensities.
     * 
     * @param inParameters reference to parameters input via GUI.
     * @param wavelength center wavelength in Angstrom.
     * @param wavelengthWidth width of the wavelegth slice to be calculated (Angstrom).
     * @return total calculated intensity within the wavelength interval.
     */
    protected abstract double getContiniumIntensity(Inparameters inParameters,
            double wavelength, double wavelengthWidth);

    // This method is implemented by the classes extending this class
    /**
     * This method which is implemented by the classes extending this class
     * calculates intensities of the characteristic lines of the x-ray tube.
     * Intensity is returned as an integrated intensity within a wavelength
     * interval. The interval width of each calculated intensity is equal to
     * the natural width of the line.
     * 
     * @param inParameters reference to parameters input via GUI.
     * @param outputData an Xraypectrum object containing the calculated values.
     */
    protected abstract void calculateTubeLineIntensities(
            Inparameters inParameters,
            XraySpectrum outputData);
    
    /**
     * Method converts a width in eV to a width in Angstrom.
     * @param lineEnergy line energy in keV.
     * @param eWidth width in eV.
     * @return width in Angstrom.
     */
    protected static double getLineWidth(double lineEnergy, double eWidth) {
        // get half of line width in keV
        double halfEwidth = 0.5d * 0.001d * eWidth;
        // Calculate line width in Angstrom
        return (halfEwidth);
    }
    
   /**
    * Method transfer a transfer factor for the tube window.
    * @param atomicNumber atomic number
    * @param wavelength wavelength in Angstrom
    * @param thickness window thickness in micrometeter
    * @return transfer factor to multiply intensity
    */
    private static double getWindowTransferFactor(
            int atomicNumber,
            double tubevoltage, //tubevoltage
            double thickness) {  // micrometer
        // Get thickness in cm
        double cmThickness = 0.0001d * thickness;//Thickness unit um
        double attCoeff = AbsCoefficient
                .getAttenuationCoefficient(atomicNumber, tubevoltage);
        return Math.exp(-attCoeff * cmThickness);
    }
    
    /**
     * Method calculates a transfer factor for a tube primary filter.
     * 
     * @param filterElements List of FilterElements describing the filter components.
     * @param tubevoltage tubevoltage
     * @param thickness filter thickness in micrometer.
     * @return transfer factor to multiply intensity
     */
    private static double getFilterTransferFactor(
            List<FilterElement> filterElements,
            double tubevoltage,
            double thickness) {  // micrometer
        // Calculate attenuation coefficient
        double attCoeff = filterElements.stream()
                .map(fElement -> {
                    int atomicNumber = fElement.getSelectedElement().getAtomicNumber();
                    double conc = fElement.getConc();
                    double attC = AbsCoefficient.getAttenuationCoefficient(atomicNumber, tubevoltage);
                    return attC * conc;
                })
                .reduce(0.0d, (a, b) -> a + b);
        double cmThickness = 0.0001d * thickness;
        return Math.exp(-attCoeff * cmThickness);
    }
    
   /**
    * Method adjusts intensities depending on tube window and tube filter attenuation.
    * @param inParameters parameters of tube window and filter input by operatior.
    * @param outputData adjusted tube spectrum intensities in XraySpectrum object.
    */
    private static void windowFilterAdjustment(
            Inparameters inParameters,
            XraySpectrum outputData) {
        // Find out if there is a filter
        boolean isFilter = true;
        List<FilterElement> filterElems = inParameters.getFilterElements();
        double filterThickness = inParameters.getFilterThickness();
        if (filterElems.isEmpty() || (filterThickness == 0)) isFilter = false;
        final boolean isFilt = isFilter;
        // Get tube window atomic number and thickness
        int windowZ = inParameters.getWindowElement().getAtomicNumber();
        double windowThickness = inParameters.getWindowThickness();
        // First adjust the tube lines
        outputData.getTubeLines().stream()
                .forEach(sPart ->{
                    double tubevoltage = sPart.getTubevoltage();
                    double wTrans = getWindowTransferFactor(windowZ, tubevoltage, windowThickness);
                    double fTrans = 1.0d;
                    if (isFilt) {
                        fTrans = getFilterTransferFactor(filterElems, tubevoltage, filterThickness);
                    }
                    sPart.setIntensity(sPart.getIntensity() * wTrans * fTrans);
                });
        // Adjust continium slices
         outputData.getContinuum().stream()
                .forEach(sPart ->{
                    double tubevoltage = sPart.getTubevoltage();
                    double wTrans = getWindowTransferFactor(windowZ, tubevoltage, windowThickness);
                    double fTrans = 1.0d;
                    if (isFilt) {
                        fTrans = getFilterTransferFactor(filterElems, tubevoltage, filterThickness);
                    }
                    sPart.setIntensity(sPart.getIntensity() * wTrans * fTrans);
                });
    }
    
    /**
     * Method normalizes calculated spectrum with a factor depending on sn
     * input parameter and the max integrated intensity of a charcteristic line.
     * 
     * @param outputData adjusted tube spectrum intensities in XraySpectrum object.
     * @param maxPeak value of max integrated peak intensity used for normalization.
     */
    private static void normalizeIntensities(XraySpectrum outputData, double maxPeak) {
        
        // Get max integrated tube line intensity
        double mPeak = outputData.getTubeLines().stream()
                .map(sPart -> sPart.getIntensity())
                .max(Double::compare)
                .get();
        
        // Get factor for normalizing
        double normFac = maxPeak / mPeak;
     
        // Normalize peaks
        outputData.getTubeLines().stream()
                .forEach(sPart -> sPart.setIntensity(normFac * sPart.getIntensity()));
        
        // Normalize continium slices
         outputData.getContinuum().stream()
                .forEach(sPart -> sPart.setIntensity(normFac * sPart.getIntensity()));
    }
}
