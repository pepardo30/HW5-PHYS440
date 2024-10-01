
import XCTest

// Test class for the 1D Schrödinger Equation Solver using Runge-Kutta 4
final class SchrodingerEquationTests: XCTestCase {

    // Test case to check if the potential is loaded correctly
    func testLoadPotential() {
        let filePath = "test_potential.txt" // Assume this file is created for testing purposes
        let potential = loadPotential(from: filePath)
        
        XCTAssertFalse(potential.isEmpty, "Potential should not be empty.")
    }

    // Test case to verify that the Schrödinger equations return valid derivatives
    func testSchrodingerEquations() {
        let psi = 1.0
        let phi = 0.5
        let Vx = 1.0
        let E = 1.0
        let result = schrodingerEquations(psi: psi, phi: phi, Vx: Vx, E: E)
        
        XCTAssertNotNil(result.dPsi, "dPsi should not be nil.")
        XCTAssertNotNil(result.dPhi, "dPhi should not be nil.")
    }

    // Test case to verify that the Runge-Kutta 4 method works correctly
    func testRungeKuttaStep() {
        let psi = 1.0
        let phi = 0.5
        let Vx = 1.0
        let E = 1.0
        let dx = 0.01
        
        let result = rungeKutta4Step(psi: psi, phi: phi, Vx: Vx, E: E, dx: dx)
        
        XCTAssertNotEqual(result.newPsi, psi, "newPsi should change after one RK4 step.")
        XCTAssertNotEqual(result.newPhi, phi, "newPhi should change after one RK4 step.")
    }
    
    // Test case to solve the Schrödinger equation for a simple potential
    func testSolveSchrodinger() {
        let potential = [1.0, 2.0, 1.5, 1.0, 0.5] // Simple test potential
        let E = 1.0
        let initialPsi = 0.0
        let initialPhi = 1.0
        let dx = 0.01
        
        let wavefunction = solveSchrodinger(potential: potential, E: E, initialPsi: initialPsi, initialPhi: initialPhi, dx: dx)
        
        XCTAssertEqual(wavefunction.count, potential.count, "The wavefunction should have the same number of points as the potential.")
        XCTAssertNotEqual(wavefunction[1], initialPsi, "The wavefunction should change from the initial value after one step.")
    }
}

// Run the tests
SchrodingerEquationTests.defaultTestSuite.run()
