import numpy as np

from cromosim.micro import projection

def test_uzawa_projection():
    xyrv = np.array(
        [
            [0.0, 0.0, 0.5, 1.0],
            [0.8, 0.0, 0.5, 1.0],
        ]
    )
    contacts = np.array([[0.0, 1.0, -0.2, 1.0, 0.0]])
    desired_velocity = np.array(
        [
            [1.0, 0.0],
            [0.0, -1.0],
        ]
    )
    mass = np.array([1.0, 1.0])

    dt = 0.1
    
    info, matrix, velocity, multipliers, pressure, _ = projection(
        dt,
        xyrv,
        contacts,
        desired_velocity,
        mass,
        dmin_people=0.0,
        nb_iter_max=1_000_000,
        rho=0.1,
        tol=1e-10,
        method="uzawa",
    )

    new_xyrv = np.array(
        [
            [0.0+dt*velocity[0,0], 0.0+dt*velocity[0,1], 0.5, 1.0],
            [0.8+dt*velocity[1,0], 0.0+dt*velocity[1,1], 0.5, 1.0],
        ]
    )

    print(f"xy={xyrv[:,:2]}")
    print(f"info={info}")
    print(f"velocity={velocity}")
    print(f"multipliers={multipliers}")
    print(f"pressure={pressure}")
    print(f"new xy={new_xyrv[:,:2]}")
        
    assert info > 0
    assert matrix.shape == (1, 4)
    # np.testing.assert_allclose(velocity, [[-1.5, 0.0], [0.5, 0.0]], atol=2e-8)
    # np.testing.assert_allclose(multipliers, [1.5], atol=2e-8)
    np.testing.assert_allclose(new_xyrv[:,:2], np.array([[-0.05,  0.],[0.95,-0.1]]), atol=1e-2)
    assert pressure.shape == (2,)



def test_nesterov_projection():
    xyrv = np.array(
        [
            [0.0, 0.0, 0.5, 1.0],
            [0.8, 0.0, 0.5, 1.0],
        ]
    )
    contacts = np.array([[0.0, 1.0, -0.2, 1.0, 0.0]])
    desired_velocity = np.array(
        [
            [1.0, 0.0],
            [0.0, -1.0],
        ]
    )
    mass = np.array([1.0, 1.0])

    dt = 0.1
    
    info, matrix, velocity, multipliers, pressure, _ = projection(
        dt,
        xyrv,
        contacts,
        desired_velocity,
        mass,
        dmin_people=0.0,
        nb_iter_max=1_000_000,
        rho=0.1,
        tol=1e-10,
        method="nesterov",
    )

    new_xyrv = np.array(
        [
            [0.0+dt*velocity[0,0], 0.0+dt*velocity[0,1], 0.5, 1.0],
            [0.8+dt*velocity[1,0], 0.0+dt*velocity[1,1], 0.5, 1.0],
        ]
    )

    print(f"xy={xyrv[:,:2]}")
    print(f"info={info}")
    print(f"velocity={velocity}")
    print(f"multipliers={multipliers}")
    print(f"pressure={pressure}")
    print(f"new xy={new_xyrv[:,:2]}")
        
    assert info > 0
    assert matrix.shape == (1, 4)
    # np.testing.assert_allclose(velocity, [[-1.5, 0.0], [0.5, 0.0]], atol=2e-8)
    # np.testing.assert_allclose(multipliers, [1.5], atol=2e-8)
    np.testing.assert_allclose(new_xyrv[:,:2], np.array([[-0.05,  0.],[0.95,-0.1]]), atol=1e-2)
    assert pressure.shape == (2,)



def test_cppnesterov_projection():
    xyrv = np.array(
        [
            [0.0, 0.0, 0.5, 1.0],
            [0.8, 0.0, 0.5, 1.0],
        ]
    )
    contacts = np.array([[0.0, 1.0, -0.2, 1.0, 0.0]])
    desired_velocity = np.array(
        [
            [1.0, 0.0],
            [0.0, -1.0],
        ]
    )
    mass = np.array([1.0, 1.0])

    dt = 0.1
    
    info, matrix, velocity, multipliers, pressure, _ = projection(
        dt,
        xyrv,
        contacts,
        desired_velocity,
        mass,
        dmin_people=0.0,
        nb_iter_max=1_000_000,
        rho=0.1,
        tol=1e-10,
        method="cppnesterov",
    )

    new_xyrv = np.array(
        [
            [0.0+dt*velocity[0,0], 0.0+dt*velocity[0,1], 0.5, 1.0],
            [0.8+dt*velocity[1,0], 0.0+dt*velocity[1,1], 0.5, 1.0],
        ]
    )

    print(f"xy={xyrv[:,:2]}")
    print(f"info={info}")
    print(f"velocity={velocity}")
    print(f"multipliers={multipliers}")
    print(f"pressure={pressure}")
    print(f"new xy={new_xyrv[:,:2]}")
        
    assert info > 0
    assert matrix.shape == (1, 4)
    # np.testing.assert_allclose(velocity, [[-1.5, 0.0], [0.5, 0.0]], atol=2e-8)
    # np.testing.assert_allclose(multipliers, [1.5], atol=2e-8)
    np.testing.assert_allclose(new_xyrv[:,:2], np.array([[-0.05,  0.],[0.95,-0.1]]), atol=1e-2)
    assert pressure.shape == (2,)

