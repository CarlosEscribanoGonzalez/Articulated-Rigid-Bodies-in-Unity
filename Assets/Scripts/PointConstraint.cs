using UnityEngine;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;

/// <summary>
/// Basic point constraint between two rigid bodies.
/// </summary>
public class PointConstraint : MonoBehaviour, IConstraint
{
    /// <summary>
    /// Default constructor. All zero. 
    /// </summary>
    public PointConstraint()
    {
        Manager = null;
    }

    #region EditorVariables

    public float Stiffness;

    public RigidBody bodyA;
    public RigidBody bodyB;

    #endregion

    #region OtherVariables

    int index;
    private PhysicsManager Manager;

    protected Vector3 pointA;
    protected Vector3 pointB;
    private MatrixXD I;
    Vector3 PosA => (bodyA != null) ? bodyA.PointLocalToGlobal(pointA) : pointA;
    Vector3 PosB => (bodyB != null) ? bodyB.PointLocalToGlobal(pointB) : pointB;
    Vector3 Xa => (bodyA != null) ? bodyA.m_pos : pointA;
    Vector3 Xb => (bodyB != null) ? bodyB.m_pos : pointB;

    #endregion

    #region MonoBehaviour

    private void Start()
    {
        I = DenseMatrixXD.CreateIdentity(3);
    }

    // Update is called once per frame
    void Update()
    {
        Vector3 pos = 0.5f * (PosA + PosB);

        // Apply the position
        Transform xform = GetComponent<Transform>();
        xform.position = pos;
    }

    #endregion

    #region IConstraint

    public void Initialize(int ind, PhysicsManager m)
    {
        index = ind;
        Manager = m;

        // Initialize local positions. We assume that the object is connected to a Sphere mesh.
        Transform xform = GetComponent<Transform>();
        if (xform == null)
        {
            System.Console.WriteLine("[ERROR] Couldn't find any transform to the constraint");
        }
        else
        {
            System.Console.WriteLine("[TRACE] Succesfully found transform connected to the constraint");
        }

        // Initialize kinematics
        Vector3 pos = xform.position;

        // Local positions on objects
        pointA = (bodyA != null) ? bodyA.PointGlobalToLocal(pos) : pos;
        pointB = (bodyB != null) ? bodyB.PointGlobalToLocal(pos) : pos;
    }

    public int GetNumConstraints()
    {
        return 3;
    }

    public void GetConstraints(VectorXD c)
    {
        c.SetSubVector(index, 3, ConstraintForce());
    }

    private VectorXD ConstraintForce()
    {
        return Utils.ToVectorXD(PosA - PosB);
    }

    public void GetConstraintJacobian(MatrixXD dcdx)
    {
        if(bodyA != null)
        {
            dcdx.SetSubMatrix(index, bodyA.index, GetPointJacobian_pos(ChosenPoint.PointA));
            dcdx.SetSubMatrix(index, bodyA.index + 3, GetPointJacobian_rot(ChosenPoint.PointA));
        }
        if(bodyB != null)
        {
            dcdx.SetSubMatrix(index, bodyB.index, GetPointJacobian_pos(ChosenPoint.PointB));
            dcdx.SetSubMatrix(index, bodyB.index + 3, GetPointJacobian_rot(ChosenPoint.PointB));
        }  
    }

    public void GetForce(VectorXD force)
    {
        if (bodyA != null)
        {
            //bodyA's force:
            VectorXD Fa = -Stiffness * GetPointJacobian_pos(ChosenPoint.PointA).Transpose() * ConstraintForce();
            force.SetSubVector(bodyA.index, 3, force.SubVector(bodyA.index, 3) + Fa);
            //bodyA's torque:
            VectorXD Ta = -Stiffness * GetPointJacobian_rot(ChosenPoint.PointA).Transpose() * ConstraintForce();
            force.SetSubVector(bodyA.index + 3, 3, force.SubVector(bodyA.index + 3, 3) + Ta);
        }
        if (bodyB != null)
        {
            //bodyB's force:
            VectorXD Fb = -Stiffness * GetPointJacobian_pos(ChosenPoint.PointB).Transpose() * ConstraintForce();
            force.SetSubVector(bodyB.index, 3, force.SubVector(bodyB.index, 3) + Fb);
            //bodyB's torque:
            VectorXD Tb = -Stiffness * GetPointJacobian_rot(ChosenPoint.PointB).Transpose() * ConstraintForce();
            force.SetSubVector(bodyB.index + 3, 3, force.SubVector(bodyB.index + 3, 3) + Tb);
        }
    }

    public void GetForceJacobian(MatrixXD dFdx, MatrixXD dFdv)
    {
        int indexA = bodyA != null ? bodyA.index : -1;
        int indexB = bodyB != null ? bodyB.index : -1;
        MatrixXD dCdxa = bodyA != null ? GetPointJacobian_pos(ChosenPoint.PointA) : new DenseMatrixXD(3, 3);
        MatrixXD dCdxb = bodyB != null ? GetPointJacobian_pos(ChosenPoint.PointB) : new DenseMatrixXD(3, 3);
        MatrixXD dCdwa = bodyA != null ? GetPointJacobian_rot(ChosenPoint.PointA) : new DenseMatrixXD(3, 3);
        MatrixXD dCdwb = bodyB != null ? GetPointJacobian_rot(ChosenPoint.PointB) : new DenseMatrixXD(3, 3);
        InsertInJacobian(dFdx, indexA, indexA, -Stiffness * dCdxa.Transpose() * dCdxa);
        InsertInJacobian(dFdx, indexB, indexB, -Stiffness * dCdxb.Transpose() * dCdxb);
        InsertInJacobian(dFdx, indexA, indexA + 3, -Stiffness * dCdxa.Transpose() * dCdwa);
        InsertInJacobian(dFdx, indexB, indexB + 3, -Stiffness * dCdxb.Transpose() * dCdwb);
        InsertInJacobian(dFdx, indexA, indexB, -Stiffness * dCdxa.Transpose() * dCdxb);
        InsertInJacobian(dFdx, indexB, indexA, -Stiffness * dCdxb.Transpose() * dCdxa);
        InsertInJacobian(dFdx, indexA, indexB + 3, -Stiffness * dCdxa.Transpose() * dCdwb);
        InsertInJacobian(dFdx, indexB, indexA + 3, -Stiffness * dCdxb.Transpose() * dCdwa);
        InsertInJacobian(dFdx, indexA + 3, indexA + 3, -Stiffness * dCdwa.Transpose() * dCdwa);
        InsertInJacobian(dFdx, indexB + 3, indexB + 3, -Stiffness * dCdwb.Transpose() * dCdwb);
        InsertInJacobian(dFdx, indexA + 3, indexA, -Stiffness * dCdwa.Transpose() * dCdxa);
        InsertInJacobian(dFdx, indexB + 3, indexB, -Stiffness * dCdwb.Transpose() * dCdxb);
        InsertInJacobian(dFdx, indexA + 3, indexB + 3, -Stiffness * dCdwa.Transpose() * dCdwb);
        InsertInJacobian(dFdx, indexB + 3, indexA + 3, -Stiffness * dCdwb.Transpose() * dCdwa);
        InsertInJacobian(dFdx, indexA + 3, indexB, -Stiffness * dCdwa.Transpose() * dCdxb);
        InsertInJacobian(dFdx, indexB + 3, indexA, -Stiffness * dCdwb.Transpose() * dCdxa);
    }

    private void InsertInJacobian(MatrixXD J, int forceIndex, int dofIndex, MatrixXD subMatrix)
    {
        if (forceIndex == -1 || dofIndex == -1) return;
        MatrixXD stored = J.SubMatrix(forceIndex, 3, dofIndex, 3);
        J.SetSubMatrix(forceIndex, dofIndex, stored + subMatrix);
    }

    private const float offset = 1e-4f;
    //Calculates the partial derivative of the constraint with respect to the chosen's point position
    private MatrixXD GetPointJacobian_pos(ChosenPoint chosenPoint)
    {
        Vector3[] vectorsToAdd = new Vector3[3] 
            { Vector3.right * offset, Vector3.up * offset, Vector3.forward * offset };
        VectorXD initConstraint = ConstraintForce();
        MatrixXD J = new DenseMatrixXD(3, 3);
        Vector3 initPosA = Xa;
        Vector3 initPosB = Xb;
        for(int i = 0; i < 3; i++)
        {
            if (chosenPoint == ChosenPoint.PointA) bodyA.m_pos += vectorsToAdd[i];
            else bodyB.m_pos += vectorsToAdd[i];
            VectorXD newConstraint = ConstraintForce();
            J.SetColumn(i, (newConstraint - initConstraint) / offset);
            if(bodyA != null) bodyA.m_pos = initPosA;
            if(bodyB != null) bodyB.m_pos = initPosB;
        }
        return J;
    }

    //Calculates the partial derivative of the constraint with respect to the chosen's point rotation
    public MatrixXD GetPointJacobian_rot(ChosenPoint chosenPoint)
    {
        Quaternion[] quaternionsToRotate = new Quaternion[3]
            { Utils.ToQuaternion(Vector3.right * offset), Utils.ToQuaternion(Vector3.up * offset),
                Utils.ToQuaternion(Vector3.forward * offset) };
        VectorXD initConstraint = ConstraintForce();
        MatrixXD J = new DenseMatrixXD(3, 3);
        Quaternion initRotA = bodyA != null ? bodyA.m_rot : Quaternion.identity;
        Quaternion initRotB = bodyB != null ? bodyB.m_rot : Quaternion.identity; 
        for (int i = 0; i < 3; i++)
        {
            if (chosenPoint == ChosenPoint.PointA) bodyA.m_rot = quaternionsToRotate[i] * bodyA.m_rot;
            else bodyB.m_rot = quaternionsToRotate[i] * bodyB.m_rot;
            VectorXD newConstraint = ConstraintForce();
            J.SetColumn(i, (newConstraint - initConstraint) / offset);
            if (bodyA != null) bodyA.m_rot = initRotA;
            if (bodyB != null) bodyB.m_rot = initRotB;
        }
        return J;
    }

    public enum ChosenPoint
    {
        PointA,
        PointB
    }
    #endregion
}
