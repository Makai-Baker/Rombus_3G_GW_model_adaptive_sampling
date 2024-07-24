import rombus.reduced_basis as reduced_basis


def test_highest_error():
    err_list = [[0, 0, 1], [0, 2], [100, 0]]
    rank, idx, err = reduced_basis._get_highest_error(err_list, [])
    assert rank == 2
    assert idx == 0
    assert err == 100
