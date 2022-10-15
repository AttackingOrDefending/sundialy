from sundialy.tools import SPA, SAMPA, SOLPOS


def round_to_3_decimals(results):
    new_results = []
    for item in results:
        if isinstance(item, tuple):
            new_item = []
            for sub_item in item:
                if isinstance(sub_item, tuple):
                    new_sub_item = []
                    for sub_sub_item in sub_item:
                        new_sub_item.append(round(sub_sub_item, 3) if type(sub_sub_item) in [int, float] else sub_sub_item)
                    new_item.append(tuple(new_sub_item))
                else:
                    new_item.append(round(sub_item, 3) if type(sub_item) in [int, float] else sub_item)
            new_item = tuple(new_item)
        else:
            new_item = round(item, 3) if type(item) in [int, float] else item
        new_results.append(new_item)
    return tuple(new_results)


def test_spa():
    results = SPA(100, 1, 5, 7, 0, 0, 0, 30, 90, -1000, 30, -20, 359, 1)
    assert round_to_3_decimals(results) == ((55.09, 54.116, 194.299, -22.759, 12.533, 286.957, 35.884), (-9.977, 6.164, 1.027, 11.304, ''), (100, 1, 7.292))


def test_sampa():
    results = SAMPA(2016, 3, 9, 1, 58, 19, 0, 10.1, 148.8, 100, 1000, 25)
    assert round_to_3_decimals(results) == ((15.082, 163.508, 15.08, 163.488, -4.377, 355.746, 349.829, 74.92), (0.268, 0.281, 0, 0.0, (1001.167, 0.0, 1062.697, 81.056, 96.013, 81.056), 'Total Solar Eclipse'))

    results = SAMPA(2020, 12, 14, 16, 19, 0, 0, -40.3, -67.9, 0, 10, 0)
    assert round_to_3_decimals(results) == ((17.116, 5.897, 17.11, 6.01, -23.267, 358.078, 262.564, 72.89), (0.271, 0.278, 0.013, 5.582, (1123.668, 62.718, 1132.799, 103.132, 58.894, 43.191), 'Partial Solar Eclipse'))

    results = SAMPA(2019, 1, 2, 3, 5, 55, 0, -23.923, -130.741, -100, 2000, -10)
    assert round_to_3_decimals(results) == ((84.636, 247.099, 126.629, 227.878, -14.482, 142.065, 235.143, -36.629), (0.271, 0.252, 0.231, 100.0, (409.514, 409.514, 61.04, 61.04, 22.755, 22.755), 'No Eclipse'))


def test_solpos():
    results = SOLPOS(1990, 3, 14, 12, 3, 4, 0, 20, 3, 1000, 15)
    assert round_to_3_decimals(results) == (1.082, 1.068, 22.561, 183.763, 67.432, 67.439, 1277.451, 1383.314, 1277.451, 360.963, 1073.612, 1.202, 1.009, 0.991, 71.014, -2.523, -9.288, 354.167)

    results = SOLPOS(2016, 8, 30, 15, 0, 0, 0, -70, -160, 780, 5, 70, 30, 8, 32, 0.05, 2)
    assert round_to_3_decimals(results) == (-1, -1.0, 98.971, 89.295, -9.0, -8.971, 0, 0, 0.0, 1099.814, 1621.077, 1.073, 0.76, 1.316, 238.685, 8.694, -0.446, 159.345)

    results = SOLPOS(2100, 5, 9, 9, 5, 0, 0, 23, 48, 10, 30, 120, 150)
    assert round_to_3_decimals(results) == (1.008, 0.01, 7.351, 221.59, 82.649, 82.649, 1329.488, 1340.506, 0, 133.95, 915.192, 1.205, 1.001, 0.999, 126.247, 17.423, 3.429, 46.416)

    results = SOLPOS(1999, 7, 22, 9, 45, 37, -5, 33.65, -84.43, 1006, 27, 135, 33.65)
    assert round_to_3_decimals(results) == (1.336, 1.327, 41.59, 97.033, 48.396, 48.41, 989.666, 1323.24, 1207.549, 347.175, 1181.11, 1.202, 1.037, 0.964, 199.233, 20.284, -6.422, 121.519)
