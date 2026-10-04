using UnityEngine;
using System.Collections.Generic;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;

/// <summary>
/// Basic physics manager capable of simulating a given ISimulable
/// implementation using diverse integration methods: explicit,
/// implicit, Verlet and semi-implicit.
/// </summary>
public class PhysicsManager : MonoBehaviour 
{
	/// <summary>
	/// Default constructor. Zero all. 
	/// </summary>
	public PhysicsManager()
	{
		Paused = true;
		TimeStep = 0.01f;
		Gravity = new Vector3 (0.0f, -9.81f, 0.0f);
		IntegrationMethod = Integration.Symplectic;
	}

	/// <summary>
	/// Integration method.
	/// </summary>
	public enum Integration
	{
		Symplectic = 1,
        Implicit = 2,
        SymplecticConstraints = 3,
    };

	#region InEditorVariables

	public bool Paused;
	public float TimeStep;
    public Vector3 Gravity;
    public List<GameObject> SimObjects;
    public List<GameObject> Constraints;
    public Integration IntegrationMethod;

    #endregion

    #region OtherVariables

    private List<ISimulable> m_objs;
    private List<IConstraint> m_constraints;
    private int m_numDoFs;
    private int m_numConstraints;

    VectorXD force;
    MatrixXD mass;
    MatrixXD mInv;
    VectorXD x;
    VectorXD v;
    MatrixXD J;
    VectorXD C;
    MatrixXD posJacobian;
    MatrixXD velJacobian;

    private bool IsImplicit => IntegrationMethod == Integration.Implicit;
    private bool IsSymplecticConstraints => IntegrationMethod == Integration.SymplecticConstraints;

    #endregion

    #region MonoBehaviour

    public void Start()
    {
        Time.fixedDeltaTime = Mathf.Clamp(TimeStep, 0.0001f, 0.05f);
        //Parse the simulable objects and initialize their state indices
        m_numDoFs = 0;
        m_objs = new List<ISimulable>(SimObjects.Count);

        foreach (GameObject obj in SimObjects)
        {
            ISimulable simobj = obj.GetComponent<ISimulable>();
            if (simobj != null)
            {
                m_objs.Add(simobj);

                // Initialize simulable object
                simobj.Initialize(m_numDoFs, this);

                // Retrieve pos and vel size
                m_numDoFs += simobj.GetNumDoFs();
            }
        }

        //Parse the constraints
        m_numConstraints = 0;
        m_constraints = new List<IConstraint>(Constraints.Count);

        foreach (GameObject obj in Constraints)
        {
            IConstraint constraint = obj.GetComponent<IConstraint>();
            if (constraint != null)
            {
                m_constraints.Add(constraint);

                // Initialize constraint
                constraint.Initialize(m_numConstraints, this);

                // Retrieve the number of constraints
                m_numConstraints += constraint.GetNumConstraints();
            }
        }

        force = new DenseVectorXD(m_numDoFs);
        mass = new DenseMatrixXD(m_numDoFs, m_numDoFs);
        mInv = new DenseMatrixXD(m_numDoFs, m_numDoFs);
        x = new DenseVectorXD(m_numDoFs);
        v = new DenseVectorXD(m_numDoFs);
        J = new DenseMatrixXD(m_numConstraints, m_numDoFs);
        C = new DenseVectorXD(m_numConstraints);
        posJacobian = new DenseMatrixXD(m_numDoFs, m_numDoFs); 
        velJacobian = new DenseMatrixXD(m_numDoFs, m_numDoFs);
    }

    public void Update()
	{
		if (Input.GetKeyUp (KeyCode.P))
			this.Paused = !this.Paused;
    }

    public void FixedUpdate()
    {
        if (this.Paused)
            return; // Not simulating

        // Select integration method
        switch (this.IntegrationMethod)
        {
            case Integration.Symplectic: this.StepSymplectic(); break;
            case Integration.Implicit: this.StepImplicit(); break;
            case Integration.SymplecticConstraints: this.StepSymplecticConstraints(); break;
            default:
                throw new System.Exception("[ERROR] Should never happen!");
        }
    }

    #endregion

    private void GetVariables()
    {
        //The rest of the matrices are overwritten in their Get() function,
        //so there is no need to call Clear()
        force.Clear();
        mass.Clear();
        mInv.Clear();
        if (IsImplicit) 
        {
            posJacobian.Clear();
            velJacobian.Clear();
        }
        foreach (ISimulable s in m_objs)
        {
            s.GetVelocity(v);
            s.GetMass(mass);
            s.GetMassInverse(mInv);
            s.GetForce(force);
            if (IsImplicit) s.GetForceJacobian(posJacobian, velJacobian);
        }
        foreach (IConstraint c in m_constraints)
        {
            c.GetForce(force);
            if (IsImplicit) c.GetForceJacobian(posJacobian, velJacobian);
            else if (IsSymplecticConstraints)
            {
                c.GetConstraints(C);
                c.GetConstraintJacobian(J);
            }
        }
    }

    private void UpdateVariables(bool incremental = true)
    {
        foreach (ISimulable s in m_objs)
        {
            if (incremental) s.AdvanceIncrementalPosition(x);
            else s.SetPosition(x);
            s.SetVelocity(v);
        }
    }

    /// <summary>
    /// Performs a simulation step using Symplectic integration.
    /// </summary>
    private void StepSymplectic()
	{
        GetVariables();
        v += TimeStep * mInv * force;
        x = TimeStep * v;
        UpdateVariables();
    }

    /// <summary>
    /// Performs a simulation step using Implicit integration.
    /// </summary>
    private void StepImplicit()
    {
        GetVariables();
        MatrixXD a = mass - TimeStep * velJacobian - TimeStep * TimeStep * posJacobian;
        VectorXD b = (mass - TimeStep * velJacobian) * v + TimeStep * force;
        v = a.Solve(b);
        x = TimeStep * v;
        UpdateVariables();
    }

    /// <summary>
    /// Performs a simulation step using Symplectic integration with constrained dynamics.
    /// The constraints are treated as implicit
    /// </summary>
    private void StepSymplecticConstraints()
    {
        GetVariables();
        VectorXD vUnconstrained = v + TimeStep * mInv * force; //Velocity not taking into account the constraint
        //After replacing the first equation in the second (v), lhs and rhs are obtained to calculate lambda
        MatrixXD lhs = J * mInv * J.Transpose(); 
        VectorXD rhs = J * vUnconstrained + 1 / TimeStep * C;
        VectorXD lambda = lhs.Solve(rhs);
        //lambda is used in the second equation to calculate v
        v = vUnconstrained - mInv * J.Transpose() * lambda;
        x = TimeStep * v;
        UpdateVariables();
    }
}
