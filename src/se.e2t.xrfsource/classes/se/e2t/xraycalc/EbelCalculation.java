/*
 * EbelCalculation.java
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

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import se.e2t.abscoeffcalculate.AbsCoefficient;
import se.e2t.xraycalc.AbsorptionEdges.AbsEdge;
import se.e2t.xraycalc.TubeLines.LineInfo;
import se.e2t.xrfsource.spectrumclasses.SpectrumPart;
import se.e2t.xrfsource.spectrumclasses.XraySpectrum;

/**
 *
 * @author Kent Ericsson, e2t AB.
 * 
 * Class calculates an x-ray tube spectrum accordinbg to an algorith published
 * by Horst Ebel working at The Technical University in Vienna.
 * Class implements a couple of abstract methods declared in the parent class
 * SourceCalculation.
 *
 * Publications:
 * 1. X-ray Tube Spectra. Horst Ebel, X-ray Spectrometry, Vol 28,pages 255-266.
 * 2. Lt, La12, Ln, Lb123456, Ly123 spectra of x-ray tubes. Horst Ebel,
 * X-ray Spectrometry 32, pages 46-51.
 */
public class EbelCalculation extends SourceCalculation {
    
    public EbelCalculation() {
        super();
    }

    /**
     * Method produces an intensity per Angstrom value for a certain wavelength.
     * @param inParameters reference to parameters input via GUI.
     * @param tubevoltage tubevoltage in Angstrom.
     * @param tubevoltageWidth width of the tubevoltage slice to be calculated (Angstrom).
     * @return total calculated intensity within wavelenth interval.
     */
    @Override
    protected double getContiniumIntensity(Inparameters inParameters,
            double Ei, double tubevoltageWidth) {

        // Get tube anode atomic number
        int z = inParameters.getAnodeElement().getAtomicNumber();
        double zD = (double) z;

        // Calculate the x exponent
        double energy0 = inParameters.getTubeVoltage();
        double energyi = Ei;
        double xExponent = 1.109d - 0.00435d * zD + 0.00175d * energy0;
        // Calculate variables of long expression including photelectric mass absorption
        double tauEj = AbsCoefficient.getTau(z, energyi);
        double sinPhi = Math.sin(inParameters.getInAngle() * Inparameters.ANGLE_CONV);
        double sinEpsilon = Math.sin(inParameters.getOutAngle() * Inparameters.ANGLE_CONV);
        double rouZ = getRouZ(energy0, energyi, z);
        double longExpression = tauEj * 2.0d * rouZ
                * (sinPhi / sinEpsilon);
        double fFactor = (1.0d - Math.exp(-longExpression)) / longExpression;
        double deltaE =tubevoltageWidth;
        double integratedIntensity = 1.35e9d * zD
                * Math.pow(((energy0 / energyi) - 1.0d), xExponent)
                * fFactor * deltaE;

        // Return a per Angstrom value
        return integratedIntensity;
    }

    /**
     * Method calculates the rouz variable as described in Ebels paper
     *
     * @param E0 = X-ray tube voltage in keV.
     * @param Ei = When method is called during calculation of continuum
     * intensity this is the wavelength of the continuum slice. When method is
     * called during calculation of tube line intensity then this is the
     * wavelength of the absorption edge associated with the line.
     * @param atomZ = atomic number of x-reay tube target material.
     * @return = value of rouz variable.
     */
    private double getRouZ(double E0, double Ei, int atomZ) {

        double zD = (double) atomZ;
        double j = 0.0135d * zD;
        double rouZm = (AtomicWeights.getRelAtomicWeight(atomZ) / zD)
                * (0.787e-5d * Math.sqrt(j) * Math.pow(E0, 1.5d)
                + 0.735e-6d * E0 * E0);
        double lnZ = Math.log(zD);
        double m = 0.1382d - (0.9211d / Math.sqrt(zD));
        double eta = Math.pow(E0, m) * (0.1904d - 0.2236d * lnZ
                + 0.1292d * lnZ * lnZ - 0.0149d * lnZ * lnZ * lnZ);
        double lnU0 = Math.log((E0 / Ei));
        double rouZ = rouZm * lnU0 * ((0.49269d - 1.0987d * eta + 0.78557d * eta * eta)
                / (0.70256d - 1.09865d * eta + 1.0046d * eta * eta + lnU0));
        return rouZ;
    }

    /**
     * Method calculates intensities of the characteristic lines of the x-ray tube.
     * spectrum. Intensity is stored as the total intensity of the line together
     * with a width which is the natural width of the line.
     * @param inParameters reference to parameters input via GUI.
     * @param outputData an XraySpectrum object containing the calculated values.
     */
    @Override
    protected void calculateTubeLineIntensities(Inparameters inParameters,
            XraySpectrum outputData) {

        // Create a temporary storage for the lines
        List<SpectrumPart> allLines = new ArrayList<>();

        // Get tube target atomic number
        int z = inParameters.getAnodeElement().getAtomicNumber();
        double zD = (double) z;

        // Calculate K line intensities according to Ebels first paper
        final Map<TubeLines.XrfLine, Map<Integer, TubeLines.LineInfo>> majorLineInfo
                = TubeLines.getMajorLineInfo();
        final double zK = 2.0d;
        final double bK = 0.35d;
        final double constK = 5.0e13d;

        double energy0 = inParameters.getTubeVoltage();

        majorLineInfo.keySet().stream()
                .filter(xrfLine -> TubeLines.isKline(xrfLine)) // Major lines includes some L and a M line
                .filter(xrfLine -> Optional.ofNullable(majorLineInfo.get(xrfLine).get(z)).isPresent())                .filter(xrfLine -> AbsorptionEdges.getEdge(xrfLine).isPresent())
                .filter(xrfLine -> FlourYield.getYield(z, AbsorptionEdges.getEdge(xrfLine).get()).isPresent())
                .filter(xrfLine -> TransProbabilities.getTransProb(z, xrfLine).isPresent())
                .forEach(xrfLine -> {
                    LineInfo lineInfo = majorLineInfo.get(xrfLine).get(z);
                    double Ei = lineInfo.getEnergy();
                    double edgeEnergy = lineInfo.getAbsorptionEdge();
                    // Calculate U0, the overvoltage ratio
                    double u0 = energy0 / edgeEnergy;
                    // Line exists if tube voltage is above absorption edge
                    if (u0 > 1) {
                        // Calculate components of stopping power factor
                        double firstParantesis = u0 * Math.log(u0) + 1.0d - u0;
                        double bigRoot = Math.sqrt(0.0135d * zD / edgeEnergy);
                        double nominator = Math.sqrt(u0) * Math.log(u0) + 2.0d * (1.0d - Math.sqrt(u0));
                        double squareBracket = 1.0d + 16.05d * bigRoot * (nominator / firstParantesis);
                        double sPowFactor = ((zK * bK) / zD) * firstParantesis * squareBracket;
                        // Calculate the f function
                        double tau = AbsCoefficient.getTau(z, Ei);
                        double sinPhi = Math.sin(inParameters.getInAngle() * Inparameters.ANGLE_CONV);
                        double sinEpsilon = Math.sin(inParameters.getOutAngle() * Inparameters.ANGLE_CONV);
                        double rouZ = getRouZ(inParameters.getTubeVoltage(),edgeEnergy, z);
                        double longExpression = tau * 2.0d * rouZ
                                * (sinPhi / sinEpsilon);
                        double fFunction = (1.0d - Math.exp(-longExpression)) / longExpression;
                        // Calculate r
                        double r = 1.0d - 0.0081517d * zD + 3.613e-5d * zD * zD
                                + 0.009583d * zD * Math.exp(-u0) + 0.001141d * energy0;
                        double omegaJK = FlourYield.getYield(z, AbsorptionEdges.getEdge(xrfLine).get()).get();
                        double pJKL = TransProbabilities.getTransProb(z, xrfLine).get();
                        double evwidth = lineInfo.getLineWidth();
                        // Calculate line width in Angstrom
                        double lineWidth = getLineWidth(lineInfo.getEnergy(), evwidth);
                        // Calculate according to formula
                        double intensity = constK * sPowFactor * r * omegaJK * pJKL * fFunction;

                        // Store calculated value
                        allLines.add(new SpectrumPart(Ei, lineWidth, intensity));
                    }
                });

        // Calculate L line intensities according to Ebels second paper
        final Map<TubeLines.XrfLine, Map<Integer, TubeLines.LineInfo>> lLineInfo
                = TubeLines.getLlineInfo();
        final double zL = 8.0d;
        final double bL = 0.25d;
        final double fCorr = -0.4814d + 0.03781 * zD - 2.413e-4d * zD * zD;

        lLineInfo.keySet().stream()
                .filter(xrfLine -> AbsorptionEdges.getEdge(xrfLine).isPresent())
                .filter(xrfLine -> Optional.ofNullable(lLineInfo.get(xrfLine).get(z)).isPresent()) 
                .filter(xrfLine -> FlourYield.getYield(z, AbsorptionEdges.getEdge(xrfLine).get()).isPresent())
                .filter(xrfLine -> TransProbabilities.getTransProb(z, xrfLine).isPresent())
                .forEach(xrfLine -> {
                    LineInfo lineInfo = lLineInfo.get(xrfLine).get(z);
                    double Ei = lineInfo.getEnergy();
                    double edgeEnergy = lineInfo.getAbsorptionEdge();
                    // Calculate U0, the overvoltage ratio
                    double u0 = energy0 / edgeEnergy;
                    // Line exists if tube voltage is above absorption edge
                    if (u0 > 1) {
                        // Calculate components of stopping power factor
                        double firstParantesis = u0 * Math.log(u0) + 1.0d - u0;
                        double bigRoot = Math.sqrt(0.0135d * zD / edgeEnergy);
                        double nominator = Math.sqrt(u0) * Math.log(u0) + 2.0d * (1.0d - Math.sqrt(u0));
                        double squareBracket = 1.0d + 16.05d * bigRoot * (nominator / firstParantesis);
                        double sPowFactor = ((zL * bL) / zD) * firstParantesis * squareBracket;
                        // Calculate the f function
                        double tau = AbsCoefficient.getTau(z, Ei);
                        double sinPhi = Math.sin(inParameters.getInAngle() * Inparameters.ANGLE_CONV);
                        double sinEpsilon = Math.sin(inParameters.getOutAngle() * Inparameters.ANGLE_CONV);
                        double rouZ = getRouZ(inParameters.getTubeVoltage(), edgeEnergy, z);
                        double longExpression = tau * 2.0d * rouZ
                                * (sinPhi / sinEpsilon);
                        double fFunction = (1.0d - Math.exp(-longExpression)) / longExpression;
                        // Calculate r
                        double r = 1.0d - 0.0081517d * zD + 3.613e-5d * zD * zD
                                + 0.009583d * zD * Math.exp(-u0) + 0.001141 * energy0;
                        double omegaJK = FlourYield.getYield(z, AbsorptionEdges.getEdge(xrfLine).get()).get();
                        double pJKL = TransProbabilities.getTransProb(z, xrfLine).get();
                        double evwidth = lineInfo.getLineWidth();
                        // Calculate line width in Angstrom
                        double lineWidth = getLineWidth(lineInfo.getEnergy(), evwidth);
                        // Calculate the Const factor
                        double constX = 4.94e13d;
                        AbsEdge absEdgeType = AbsorptionEdges.getEdge(xrfLine).get();
                        switch (absEdgeType) {
                            case L1_EDGE:
                                constX = fCorr * 0.71e13d;
                                break;
                            case L2_EDGE:
                                constX = fCorr * 2.70e13d;
                        }
                        // Calculate line intensity
                        double intensity = constX * sPowFactor * r * omegaJK * pJKL * fFunction;
                        
                        // Store calculated intensity
                        allLines.add(new SpectrumPart(Ei, lineWidth, intensity));
                    }
                });
        // Sort lines in wavelength order and add to output data
        allLines.stream()
                .sorted(Comparator.comparing(SpectrumPart::getTubevoltage))
                .forEach(specPart -> outputData.getTubeLines().add(specPart));
    }
}
