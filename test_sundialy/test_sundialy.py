from sundialy import AnalemmaticHorizontal


def round_to_3_decimals(results):
    if type(results) in [int, float]:
        return round(results, 3)
    new_results = {}
    for key, value in results.items():
        if isinstance(value, tuple):
            new_value = []
            for sub_value in value:
                new_value.append(round(sub_value, 3) if type(sub_value) in [int, float] else sub_value)
            new_results[key] = tuple(new_value)
        else:
            new_results[key] = round(value, 3) if type(value) in [int, float] else value
    return new_results


def test_sundial():
    sundial = AnalemmaticHorizontal(latitude=34, longitude=-118, timezone=-7, show=False)
    sundial.create_sundial()
    assert round_to_3_decimals(sundial.height) == 2.796
    assert round_to_3_decimals(sundial.gnomon_movement) == {'Jan 1': -0.88, 'Feb 1': -0.637, 'Mar 1': -0.275, 'Apr 1': 0.165, 'May 1': 0.559, 'Jun 1': 0.84, 'Jun 21': 0.898, 'Jul 1': 0.884, 'Aug 1': 0.673, 'Sep 1': 0.301, 'Oct 1': -0.116, 'Nov 1': -0.534, 'Dec 1': -0.829, 'Dec 21': -0.898}
    assert round_to_3_decimals(sundial.hour_locations) == {0: (-180.0, 0.0, -1.398), 1: (-154.398, -0.646, -1.349), 2: (-134.085, -1.249, -1.209), 3: (-119.214, -1.768, -0.989), 4: (-107.893, -2.163, -0.698), 5: (-98.522, -2.412, -0.361), 6: (-90.0, -2.5, 0.0), 7: (-81.478, -2.412, 0.361), 8: (-72.107, -2.163, 0.698), 9: (-60.786, -1.768, 0.989), 10: (-45.915, -1.249, 1.209), 11: (-25.602, -0.646, 1.349), 12: (0.0, 0.0, 1.398), 13: (25.602, 0.646, 1.349), 14: (45.915, 1.249, 1.209), 15: (60.786, 1.768, 0.989), 16: (72.107, 2.163, 0.698), 17: (81.478, 2.412, 0.361), 18: (90.0, 2.5, 0.0), 19: (98.522, 2.412, -0.361), 20: (107.893, 2.163, -0.698), 21: (119.214, 1.768, -0.989), 22: (134.085, 1.249, -1.209), 23: (154.398, 0.646, -1.349)}
    assert round_to_3_decimals(sundial.significant_eot) == {'Jan 1': -55.298, 'Feb 1': -65.469, 'Feb 11': -66.189, 'Mar 1': -64.399, 'Apr 1': -56.012, 'May 1': -49.168, 'May 14': -48.361, 'Jun 1': -49.788, 'Jul 1': -55.812, 'Jul 26': -58.563, 'Aug 1': -58.403, 'Sep 1': -52.182, 'Oct 1': -41.831, 'Nov 1': -35.595, 'Nov 3': -35.557, 'Dec 1': -40.813}

    sundial = AnalemmaticHorizontal(latitude=-34, longitude=118, correct_for_longitude=True, timezone=7, years=[2022, 2023, 2024], show=False)
    sundial.create_sundial()
    assert round_to_3_decimals(sundial.height) == 2.796
    assert round_to_3_decimals(sundial.gnomon_movement) == {'Jan 1': 0.882, 'Feb 1': 0.646, 'Mar 1': 0.287, 'Apr 1': -0.158, 'May 1': -0.553, 'Jun 1': -0.837, 'Jun 21': -0.898, 'Jul 1': -0.885, 'Aug 1': -0.678, 'Sep 1': -0.307, 'Oct 1': 0.109, 'Nov 1': 0.528, 'Dec 1': 0.826, 'Dec 21': 0.898}
    assert round_to_3_decimals(sundial.hour_locations) == {0: (157.566, 0.562, -1.361), 1: (136.443, 1.174, -1.234), 2: (120.949, 1.703, -1.021), 3: (109.261, 2.118, -0.74), 4: (99.702, 2.391, -0.409), 5: (91.119, 2.496, -0.049), 6: (82.644, 2.433, 0.314), 7: (73.441, 2.207, 0.656), 8: (62.46, 1.826, 0.952), 9: (48.175, 1.323, 1.184), 10: (28.667, 0.731, 1.337), 11: (3.573, 0.087, 1.396), 12: (-22.434, -0.562, 1.361), 13: (-43.557, -1.174, 1.234), 14: (-59.051, -1.703, 1.021), 15: (-70.739, -2.118, 0.74), 16: (-80.298, -2.391, 0.409), 17: (-88.881, -2.496, 0.049), 18: (-97.356, -2.433, -0.314), 19: (-106.559, -2.207, -0.656), 20: (-117.54, -1.826, -0.952), 21: (-131.825, -1.323, -1.184), 22: (-151.333, -0.731, -1.337), 23: (-176.427, -0.087, -1.396)}
    assert round_to_3_decimals(sundial.significant_eot) == {'Jan 1': -3.187, 'Feb 1': -13.433, 'Feb 11': -14.182, 'Mar 1': -12.378, 'Apr 1': -3.977, 'May 1': 2.852, 'May 14': 3.65, 'Jun 1': 2.208, 'Jul 1': -3.82, 'Jul 26': -6.56, 'Aug 1': -6.386, 'Sep 1': -0.141, 'Oct 1': 10.211, 'Nov 1': 16.417, 'Nov 3': 16.453, 'Dec 1': 11.157}

    sundial = AnalemmaticHorizontal(latitude=0, longitude=0, correct_for_longitude=True, show=False)
    sundial.create_sundial("sundial.jpg", "corrections.jpg")
    assert round_to_3_decimals(sundial.height) == 0.
    assert round_to_3_decimals(sundial.gnomon_movement) == {'Jan 1': -1.062, 'Feb 1': -0.772, 'Mar 1': -0.337, 'Apr 1': 0.194, 'May 1': 0.67, 'Jun 1': 1.011, 'Jun 21': 1.084, 'Jul 1': 1.067, 'Aug 1': 0.816, 'Sep 1': 0.367, 'Oct 1': -0.135, 'Nov 1': -0.64, 'Dec 1': -0.998, 'Dec 21': -1.084}
    assert round_to_3_decimals(sundial.hour_locations) == {0: (-180.0, 0.0, -0.0), 1: (-90.0, -0.646, -0.0), 2: (-90.0, -1.249, -0.0), 3: (-90.0, -1.768, -0.0), 4: (-90.0, -2.163, -0.0), 5: (-90.0, -2.412, -0.0), 6: (-90.0, -2.5, 0.0), 7: (-90.0, -2.412, 0.0), 8: (-90.0, -2.163, 0.0), 9: (-90.0, -1.768, 0.0), 10: (-90.0, -1.249, 0.0), 11: (-90.0, -0.646, 0.0), 12: (0.0, 0.0, 0.0), 13: (90.0, 0.646, 0.0), 14: (90.0, 1.249, 0.0), 15: (90.0, 1.768, 0.0), 16: (90.0, 2.163, 0.0), 17: (90.0, 2.412, 0.0), 18: (90.0, 2.5, 0.0), 19: (90.0, 2.412, -0.0), 20: (90.0, 2.163, -0.0), 21: (90.0, 1.768, -0.0), 22: (90.0, 1.249, -0.0), 23: (90.0, 0.646, -0.0)}
    assert round_to_3_decimals(sundial.significant_eot) == {'Jan 1': -3.298, 'Feb 1': -13.469, 'Feb 11': -14.189, 'Mar 1': -12.399, 'Apr 1': -4.012, 'May 1': 2.832, 'May 14': 3.639, 'Jun 1': 2.212, 'Jul 1': -3.812, 'Jul 26': -6.563, 'Aug 1': -6.403, 'Sep 1': -0.182, 'Oct 1': 10.169, 'Nov 1': 16.405, 'Nov 3': 16.443, 'Dec 1': 11.187}
    assert str(sundial) == "Width: 5, Height: 0.0\nLongitude correction: 0\nFor Jan 1, the gnomon should move -1.062 meters forwards\nFor Feb 1, the gnomon should move -0.772 meters forwards\nFor Mar 1, the gnomon should move -0.337 meters forwards\nFor Apr 1, the gnomon should move 0.194 meters forwards\nFor May 1, the gnomon should move 0.67 meters forwards\nFor Jun 1, the gnomon should move 1.011 meters forwards\nFor Jun 21, the gnomon should move 1.084 meters forwards\nFor Jul 1, the gnomon should move 1.067 meters forwards\nFor Aug 1, the gnomon should move 0.816 meters forwards\nFor Sep 1, the gnomon should move 0.367 meters forwards\nFor Oct 1, the gnomon should move -0.135 meters forwards\nFor Nov 1, the gnomon should move -0.64 meters forwards\nFor Dec 1, the gnomon should move -0.998 meters forwards\nFor Dec 21, the gnomon should move -1.084 meters forwards\nThe angle for the 0th hour is -180.0 and the coordinates: (0.0, -0.0)\nThe angle for the 1th hour is -90.0 and the coordinates: (-0.646, -0.0)\nThe angle for the 2th hour is -90.0 and the coordinates: (-1.249, -0.0)\nThe angle for the 3th hour is -90.0 and the coordinates: (-1.768, -0.0)\nThe angle for the 4th hour is -90.0 and the coordinates: (-2.163, -0.0)\nThe angle for the 5th hour is -90.0 and the coordinates: (-2.412, -0.0)\nThe angle for the 6th hour is -90.0 and the coordinates: (-2.5, 0.0)\nThe angle for the 7th hour is -90.0 and the coordinates: (-2.412, 0.0)\nThe angle for the 8th hour is -90.0 and the coordinates: (-2.163, 0.0)\nThe angle for the 9th hour is -90.0 and the coordinates: (-1.768, 0.0)\nThe angle for the 10th hour is -90.0 and the coordinates: (-1.249, 0.0)\nThe angle for the 11th hour is -90.0 and the coordinates: (-0.646, 0.0)\nThe angle for the 12th hour is 0.0 and the coordinates: (0.0, 0.0)\nThe angle for the 13th hour is 90.0 and the coordinates: (0.646, 0.0)\nThe angle for the 14th hour is 90.0 and the coordinates: (1.249, 0.0)\nThe angle for the 15th hour is 90.0 and the coordinates: (1.768, 0.0)\nThe angle for the 16th hour is 90.0 and the coordinates: (2.163, 0.0)\nThe angle for the 17th hour is 90.0 and the coordinates: (2.412, 0.0)\nThe angle for the 18th hour is 90.0 and the coordinates: (2.5, 0.0)\nThe angle for the 19th hour is 90.0 and the coordinates: (2.412, -0.0)\nThe angle for the 20th hour is 90.0 and the coordinates: (2.163, -0.0)\nThe angle for the 21th hour is 90.0 and the coordinates: (1.768, -0.0)\nThe angle for the 22th hour is 90.0 and the coordinates: (1.249, -0.0)\nThe angle for the 23th hour is 90.0 and the coordinates: (0.646, -0.0)\nFor Jan 1, the equation of time (EOT) is -3.298 minutes\nFor Feb 1, the equation of time (EOT) is -13.469 minutes\nFor Feb 11, the equation of time (EOT) is -14.189 minutes\nFor Mar 1, the equation of time (EOT) is -12.399 minutes\nFor Apr 1, the equation of time (EOT) is -4.012 minutes\nFor May 1, the equation of time (EOT) is 2.832 minutes\nFor May 14, the equation of time (EOT) is 3.639 minutes\nFor Jun 1, the equation of time (EOT) is 2.212 minutes\nFor Jul 1, the equation of time (EOT) is -3.812 minutes\nFor Jul 26, the equation of time (EOT) is -6.563 minutes\nFor Aug 1, the equation of time (EOT) is -6.403 minutes\nFor Sep 1, the equation of time (EOT) is -0.182 minutes\nFor Oct 1, the equation of time (EOT) is 10.169 minutes\nFor Nov 1, the equation of time (EOT) is 16.405 minutes\nFor Nov 3, the equation of time (EOT) is 16.443 minutes\nFor Dec 1, the equation of time (EOT) is 11.187 minutes"
