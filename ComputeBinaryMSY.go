package main

import (
	"math"
	"fmt"

	"github.com/keyan/goent/discrete"

    "log"
    "os"

    "strconv"
)

func alpha(w2 float64, w1 float64, a1 float64, phi float64, psi float64, chi float64) float64{
	nom := math.Exp(phi*w2*w1 + psi*w2*a1 + chi*w2*w1*a1)

	denom := 0.0
	RandVarValue := [2]float64{-1, 1}
	for w3 := 0; w3 < 2; w3++ {
		denom += math.Exp(phi*RandVarValue[w3]*w1 + psi*RandVarValue[w3]*a1 + chi*RandVarValue[w3]*w1*a1)
	}

	return nom / denom
}

func beta(s1 float64, w1 float64, zeta float64) float64{
	nom := math.Exp(zeta*s1*w1)

	denom := 0.0
	RandVarValue := [2]float64{-1, 1}
	for s3 := 0; s3 < 2; s3++ {
		denom += math.Exp(zeta*RandVarValue[s3]*w1)
	}

	return nom / denom
}

func pi(a1 float64, s1 float64, mu float64) float64{
	nom := math.Exp(mu*a1*s1)

	denom := 0.0
	RandVarValue := [2]float64{-1, 1}
	for a2 := 0; a2 < 2; a2++ {
		denom += math.Exp(mu*RandVarValue[a2]*s1)
	}

	return nom / denom
}

func ro(w1 float64, tau float64) float64{
	nom := math.Exp(tau*w1)

	denom := 0.0
	RandVarValue := [2]float64{-1, 1}
	for w2 := 0; w2 < 2; w2++ {
		denom += math.Exp(tau*RandVarValue[w2])
	}

	return nom / denom
}

func normalise3D(p [][][]float64) [][][]float64 {
	sum := 0.0
	for i := 0; i < len(p); i++ {
		for j := 0; j < len(p[i]); j++ {
			for k := 0; k < len(p[i][j]); k++ {
				sum += p[i][j][k]
			}
		}
	}

	for i := 0; i < len(p); i++ {
		for j := 0; j < len(p[i]); j++ {
			for k := 0; k < len(p[i][j]); k++ {
				p[i][j][k] /= sum
			}
		}
	}
	return p
}

func MorphologicalComputationSYBinary(phi float64, psi float64, chi float64, zeta float64, mu float64, tau float64) float64{
	iterations := 100
	pw2w1a1 := discrete.Create3D(2, 2, 2)
	RandVarValue := [2]float64{-1, 1}
	for w2 := 0; w2 < 2; w2+=1 {
		for w1 := 0; w1 < 2; w1+=1 {
			for a1 := 0; a1 < 2; a1+=1 {
				pw2w1a1[w2][w1][a1] = 0
				for s1 := 0; s1 < 2; s1+=1 {
					pw2w1a1[w2][w1][a1] += alpha(RandVarValue[w2], RandVarValue[w1], RandVarValue[a1], phi, psi, chi) * 
											beta(RandVarValue[s1], RandVarValue[w1], zeta) * 
											pi(RandVarValue[a1], RandVarValue[s1], mu) * 
											ro(RandVarValue[w1], tau)
				}
			}
		}
	}
	pw2w1a1 = normalise3D(pw2w1a1)
	return discrete.MorphologicalComputationSY(pw2w1a1, iterations, false)
}

func GenerateMCSYBinaryOutput(phiRange []float64, psiRange []float64, chi float64) [][]float64{
	output := discrete.Create2D(len(phiRange), len(psiRange))
	for i := 0; i < len(phiRange); i++ {
		for j := 0; j < len(psiRange); j++ {
			output[i][j] = MorphologicalComputationSYBinary(phiRange[i], psiRange[j], chi, 0, 0, 0)
		}
	}
	return output
}


func ArrayOfFloatsToString(ArrayOfFloats []float64) string{
	str := ""
	for _, num := range ArrayOfFloats{
		s := (strconv.FormatFloat(num, 'E', 3, 64))
		str += s + " "
	}
	return str + "\n"
}

func write2DArray(data [][]float64, filename string) {
	file, err := os.Create("outputs/"+ filename + ".txt")
    if err != nil {
        log.Fatal("Cannot create file", err)
    }
    defer file.Close()

    outputString := ""
    for _, line := range data {
    	outputString += ArrayOfFloatsToString(line)
    }

    fmt.Fprintf(file, outputString)
}

func generateSequence(min float64, max float64, increment float64) []float64{
	length := int(float64(max - min) / increment)
	output := make([]float64, length, length)
	for i := 0; i < length; i++ {
		output[i] = min + increment*float64(i)
	}
	return output
}

func main() {
	var increment float64
	var minPhi, maxPhi, minPsi, maxPsi float64
	var chi float64
	var filename string

	if len(os.Args) != 8 {
		log.Fatal("Incorrect Number of Arguments: Expected 8 but got %d.", len(os.Args))
	}

	minPhi, err1 := strconv.ParseFloat(os.Args[1], 64)
	maxPhi, err2 := strconv.ParseFloat(os.Args[2], 64)
	minPsi, err3 := strconv.ParseFloat(os.Args[3], 64)
	maxPsi, err4 := strconv.ParseFloat(os.Args[4], 64)
	chi, err5 := strconv.ParseFloat(os.Args[5], 64) 
	increment, err6 := strconv.ParseFloat(os.Args[6], 64)
	filename = os.Args[7]


	errs := [6]error{err1, err2, err3, err4, err5, err6}
	for i := 0; i < 6; i++ {
		if errs[i] != nil {
			log.Fatal("Something is wrong with input format.")
		}
	}

	MCSYBinaryOutput := GenerateMCSYBinaryOutput(generateSequence(minPhi, maxPhi, increment),
													generateSequence(minPsi, maxPsi, increment), chi)
	write2DArray(MCSYBinaryOutput, filename)
	
}