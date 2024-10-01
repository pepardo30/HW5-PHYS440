
import Foundation

// Constants
let hbar = 1.0 // Reduced Planck constant (can be set to 1 in natural units)
let mass = 1.0 // Particle mass (can be set to 1 in natural units)
let stepSize = 0.01 // Step size for RK4
let maxSteps = 1000 // Maximum number of steps

// Function to read the potential from a file
func loadPotential(from filePath: String) -> [Double] {
    do {
        let data = try String(contentsOfFile: filePath, encoding: .utf8)
        return data.split(separator: "
").compactMap { Double($0) }
    } catch {
        print("Error loading file: \(error)")
        return []
    }
}

// Schrödinger equation system of ODEs
// dψ/dx = φ
// dφ/dx = (2m/ħ²)(V(x) - E)ψ
func schrodingerEquations(psi: Double, phi: Double, Vx: Double, E: Double) -> (dPsi: Double, dPhi: Double) {
    let dPsi = phi
    let dPhi = 2 * mass * (Vx - E) * psi / (hbar * hbar)
    return (dPsi, dPhi)
}

// Runge-Kutta 4th Order method for solving the differential equation
func rungeKutta4Step(psi: Double, phi: Double, Vx: Double, E: Double, dx: Double) -> (newPsi: Double, newPhi: Double) {
    let k1 = schrodingerEquations(psi: psi, phi: phi, Vx: Vx, E: E)
    let k2 = schrodingerEquations(psi: psi + 0.5 * k1.dPsi * dx, phi: phi + 0.5 * k1.dPhi * dx, Vx: Vx, E: E)
    let k3 = schrodingerEquations(psi: psi + 0.5 * k2.dPsi * dx, phi: phi + 0.5 * k2.dPhi * dx, Vx: Vx, E: E)
    let k4 = schrodingerEquations(psi: psi + k3.dPsi * dx, phi: phi + k3.dPhi * dx, Vx: Vx, E: E)
    
    let newPsi = psi + (dx / 6.0) * (k1.dPsi + 2 * k2.dPsi + 2 * k3.dPsi + k4.dPsi)
    let newPhi = phi + (dx / 6.0) * (k1.dPhi + 2 * k2.dPhi + 2 * k3.dPhi + k4.dPhi)
    
    return (newPsi, newPhi)
}

// Function to solve the Schrödinger equation for a given potential
func solveSchrodinger(potential: [Double], E: Double, initialPsi: Double, initialPhi: Double, dx: Double) -> [Double] {
    var psi = initialPsi
    var phi = initialPhi
    var wavefunction = [psi]
    
    for i in 0..<potential.count - 1 {
        let Vx = potential[i]
        let (newPsi, newPhi) = rungeKutta4Step(psi: psi, phi: phi, Vx: Vx, E: E, dx: dx)
        psi = newPsi
        phi = newPhi
        wavefunction.append(psi)
    }
    
    return wavefunction
}

// Example usage: Load a potential, solve the Schrödinger equation, and print the result
func main() {
    let filePath = "potential.txt" // Specify the path to the potential file
    let potential = loadPotential(from: filePath)
    
    if potential.isEmpty {
        print("No potential data loaded.")
        return
    }
    
    let energy = 1.0 // Specify the energy level to solve for
    let initialPsi = 0.0 // Initial condition for psi
    let initialPhi = 1.0 // Initial condition for phi (derivative of psi)
    
    let wavefunction = solveSchrodinger(potential: potential, E: energy, initialPsi: initialPsi, initialPhi: initialPhi, dx: stepSize)
    
    // Print or process the wavefunction
    for (x, psi) in wavefunction.enumerated() {
        print("x = \(Double(x) * stepSize), ψ(x) = \(psi)")
    }
}

main()
