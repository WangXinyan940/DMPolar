from dmpolar.wrappers.gau import check_gaussian_installation


def test_gaussian_check():
    ret = check_gaussian_installation()
    print(ret)