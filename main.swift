
import Foundation

// Main function to run the Schrödinger equation simulation
func main() {
    print("Running Schrödinger Equation Simulation...")
    
    // Example usage: Load a potential, solve the Schrödinger equation, and print the result
    let filePath = "potential.txt" // Path to the potential file
    let potential = loadPotential(from: filePath)
    
    if potential.isEmpty {
        print("No potential data loaded.")
        return
    }
    
    let energy = 1.0 // Specify the energy level to solve for
    let initialPsi = 0.0 // Initial condition for psi
    let initialPhi = 1.0 // Initial condition for phi (derivative of psi)
    let stepSize = 0.01
    
    let wavefunction = solveSchrodinger(potential: potential, E: energy, initialPsi: initialPsi, initialPhi: initialPhi, dx: stepSize)
    
    // Print or process the wavefunction
    for (x, psi) in wavefunction.enumerated() {
        print("x = \(Double(x) * stepSize), ψ(x) = \(psi)")
    }
}

main()
