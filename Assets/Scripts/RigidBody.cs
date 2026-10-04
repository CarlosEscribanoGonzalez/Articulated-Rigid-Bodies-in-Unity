using UnityEngine;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;

/// <summary>
/// Basic rigid-body model component which can be dropped onto
/// a game object.
/// </summary>
public class RigidBody : MonoBehaviour, ISimulable
{
    /// <summary>
    /// Default constructor. All zero. 
    /// </summary>
    public RigidBody()
    {
        Manager = null;
    }

    #region EditorVariables

    public float Mass;
    public float Damping;

    #endregion

    #region OtherVariables

    private PhysicsManager Manager;

    public int index;

    public Vector3 m_pos;
    public Quaternion m_rot;
    public Vector3 m_vel;
    public Vector3 m_omega;

    protected MatrixXD m_inertia0, m_inertia0inv;
    protected MatrixXD m_inertia, m_inertiainv;
    private MatrixXD I;

    #endregion

    #region MonoBehaviour

    // Update is called once per frame
    void Update()
    {
        // Apply the position and rotation to the mesh
        Transform xform = GetComponent<Transform>();
        xform.position = m_pos;
        xform.rotation = m_rot;
    }

    #endregion

    #region ISimulable

    public void Initialize(int ind, PhysicsManager m)
    {
        Manager = m;

        index = ind;

        // Initialize inertia. We assume that the object is connected to a Cube mesh.
        Transform xform = GetComponent<Transform>();
        if (xform == null)
        {
            System.Console.WriteLine("[ERROR] Couldn't find any transform to the rigid body");
        }
        else
        {
            System.Console.WriteLine("[TRACE] Succesfully found transform connected to the rigid body");
        }

        if (xform != null)
        {
            m_inertia0 = DenseMatrixXD.CreateIdentity(3);
            double[] vals;
            vals = new double[3];
            vals[0] = 1.0f / 12.0f * Mass * (xform.localScale.y * xform.localScale.y + xform.localScale.z * xform.localScale.z);
            vals[1] = 1.0f / 12.0f * Mass * (xform.localScale.x * xform.localScale.x + xform.localScale.z * xform.localScale.z);
            vals[2] = 1.0f / 12.0f * Mass * (xform.localScale.x * xform.localScale.x + xform.localScale.y * xform.localScale.y);
            m_inertia0.SetDiagonal(vals);
        }

        // Initialize kinematics
        m_pos = xform.position;
        m_rot = xform.rotation;
        m_vel = Vector3.zero;
        m_omega = Vector3.zero;

        // Initialize all inertia terms
        m_inertia0inv = m_inertia0.Inverse();
        m_inertia = Utils.WarpMatrix(m_rot, m_inertia0);
        m_inertiainv = Utils.WarpMatrix(m_rot, m_inertia0inv);
        I = DenseMatrixXD.CreateIdentity(3);
    }

    public int GetNumDoFs()
    {
        return 6;
    }

    public void GetPosition(VectorXD position)
    {
        position.SetSubVector(index, 3, Utils.ToVectorXD(m_pos));

        //Transform the quaternion to axis-angle
        Vector3 axisangle = Utils.ToAxisAngle(m_rot);
        position.SetSubVector(index + 3, 3, Utils.ToVectorXD(axisangle));
    }

    public void SetPosition(VectorXD position)
    {
        m_pos = Utils.ToVector3(position.SubVector(index, 3));

        //Transform the axis-angle to quaternion
        Vector3 axisangle = Utils.ToVector3(position.SubVector(index + 3, 3));
        m_rot = Utils.ToQuaternion(axisangle);

        m_inertia = Utils.WarpMatrix(m_rot, m_inertia0);
        m_inertiainv = Utils.WarpMatrix(m_rot, m_inertia0inv);
    }

    public void AdvanceIncrementalPosition(VectorXD position)
    {
        m_pos += Utils.ToVector3(position.SubVector(index, 3));

        //Transform the axis-angle to quaternion
        Vector3 axisangle = Utils.ToVector3(position.SubVector(index + 3, 3));
        m_rot = Utils.ToQuaternion(axisangle) * m_rot;

        m_inertia = Utils.WarpMatrix(m_rot, m_inertia0);
        m_inertiainv = Utils.WarpMatrix(m_rot, m_inertia0inv);
    }

    public void GetVelocity(VectorXD velocity)
    {
        velocity.SetSubVector(index, 3, Utils.ToVectorXD(m_vel));
        velocity.SetSubVector(index + 3, 3, Utils.ToVectorXD(m_omega));
    }

    public void SetVelocity(VectorXD velocity)
    {
        m_vel = Utils.ToVector3(velocity.SubVector(index, 3));
        m_omega = Utils.ToVector3(velocity.SubVector(index + 3, 3));
    }

    public void GetForce(VectorXD force)
    {
        Vector3 rbForce = Mass * Manager.Gravity; //Gravity force
        rbForce += -Damping * Mass * m_vel; //Linear velocity damping
        force.SetSubVector(index, 3, force.SubVector(index, 3) + Utils.ToVectorXD(rbForce));
        VectorXD rbTorque = -Damping * m_inertia * Utils.ToVectorXD(m_omega); //Angular damping
        rbTorque -= Utils.Skew(m_omega) * (m_inertia * Utils.ToVectorXD(m_omega)); //Gyroscopic term
        force.SetSubVector(index + 3, 3, force.SubVector(index + 3, 3) + rbTorque);
    }

    public void GetForceJacobian(MatrixXD dFdx, MatrixXD dFdv)
    {
        MatrixXD linearTerm = -Damping * Mass * I; 
        dFdv.SetSubMatrix(index, index, dFdv.SubMatrix(index, 3, index, 3) + linearTerm);
        MatrixXD angularTerm = -Damping * m_inertia * I;
        //Gyroscopic term has to be considered (w x (Mw)):
        MatrixXD gyroTerm = Utils.Skew(Utils.ToVector3(m_inertia * Utils.ToVectorXD(m_omega)));
        gyroTerm += Utils.Skew(m_omega) * m_inertia;
        angularTerm -= gyroTerm;
        dFdv.SetSubMatrix(index + 3, index + 3, 
            dFdv.SubMatrix(index + 3, 3, index + 3, 3) + angularTerm);
    }

    public void GetMass(MatrixXD mass)
    {
        mass.SetSubMatrix(index, index, mass.SubMatrix(index, 3, index, 3)
            + Mass * DenseMatrixXD.CreateIdentity(3));
        mass.SetSubMatrix(index + 3, index + 3, mass.SubMatrix(index + 3, 3, index + 3, 3)
            + m_inertia);
    }

    public void GetMassInverse(MatrixXD massInv)
    {
        massInv.SetSubMatrix(index, index, massInv.SubMatrix(index, 3, index, 3)
            + 1.0f / Mass * DenseMatrixXD.CreateIdentity(3));
        massInv.SetSubMatrix(index + 3, index + 3, massInv.SubMatrix(index + 3, 3, index + 3, 3)
            + m_inertiainv);
    }

    public void FixVector(VectorXD v)
    {
    }

    public void FixMatrix(MatrixXD M)
    {
    }

    #endregion

    #region OtherMethods

    public Vector3 PointGlobalToLocal(Vector3 p)
    {
        return Quaternion.Inverse(this.m_rot) * (p - this.m_pos);
    }

    public Vector3 PointLocalToGlobal(Vector3 p)
    {
        return this.m_pos + this.m_rot * p;
    }

    public Vector3 VecLocalToGlobal(Vector3 v)
    {
        return this.m_rot * v;
    }

    #endregion
}
