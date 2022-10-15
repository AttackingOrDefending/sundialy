from typing import Union, Iterable, Tuple, Optional, Dict, List
import datetime
from .tools import SPA
from shapely import affinity, geometry  # type: ignore
from matplotlib.patches import Ellipse  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import math

NUMBER_TYPE = Union[int, float]


class Horizontal:
    """
    The class that creates the analemmatic sundial.

    :param latitude: The latitude of the sundial.
    :param longitude: The longitude of the sundial.
    :param width: The width of the sundial.
    :param timezone: The timezone that the sundial is in.
    :param correct_for_longitude: Whether to add the longitude correction.
    :param years: The year(s) to calculate the data for.
    :param elevation: The elevation of the sundial.
    :param show: Whether to show the images after creating them.

    :data:`height` The height of the analemmatic sundial.

    :data:`gnomon_movement` Where the gnomon should be each day.

    :data:`hour_locations` Where each hour should be placed on the ellipse.

    :data:`significant_eot` The important Equation of Time values.

    :data:`average_eot` The Equation of Time every day.
    """

    YEAR_NOW = datetime.datetime.now().year
    MONTH_TO_NUMBER = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6,
                       'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}
    NUMBER_TO_MONTH = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                       7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}

    def __init__(self, latitude: NUMBER_TYPE, longitude: NUMBER_TYPE, width: NUMBER_TYPE = 5, timezone: NUMBER_TYPE = 0,
                 correct_for_longitude: bool = False, years: Union[int, Iterable[int]] = YEAR_NOW,
                 elevation: NUMBER_TYPE = 0, show: bool = True) -> None:
        if not show:
            plt.switch_backend("Agg")
        self.latitude = latitude if latitude else 1e-9
        self.longitude = longitude
        self.width = width
        self.height = math.sin(math.radians(abs(self.latitude))) * self.width
        self.timezone = timezone
        self.correct_for_longitude = correct_for_longitude
        self.years: Iterable[int]
        if isinstance(years, int):
            self.years = [years]
        else:
            self.years = years
        self.elevation = elevation
        self.show = show

        self.gnomon_movement: Dict[str, float] = {}
        self.hour_locations: Dict[int, Tuple[float, float, float]] = {}
        self.significant_eot: Dict[str, float] = {}
        self.average_eot: List[float] = []

    def _hour_location(self, arc: NUMBER_TYPE) -> Tuple[float, float]:
        """
        Finds the coordinates for this arc.

        :param arc: The angle with (0, x) being 0 degrees, which is then converted to the angle with (x, 0) being 0
            degrees.
        :return: The x and y coordinates of the location of the mark.
        """
        ell = geometry.Point((0, 0)).buffer(1)
        ell = affinity.scale(ell, self.width / 2, self.height / 2)
        right_angle = 180 - ((arc + 90) % 360)
        line = geometry.LineString([(0., 0.), (self.width / 2, 0)])
        line = affinity.rotate(line, right_angle, (0, 0))
        intersections = ell.intersection(line)
        coors = intersections.coords.xy
        x_ans, y_ans = coors[0][1], coors[1][1]
        return x_ans, y_ans

    def _move(self, month_number: int, day: int) -> float:
        """
        Calculates where the gnomon should move.
        Y = Width / 2 * cos(latitude) * tan(declination)

        :param month_number: The month of the year.
        :param day: The day of the month.
        :return: How many meters forward the gnomon should move.
        """
        declinations = []
        time_with_offset = datetime.datetime(2, month_number, day, 0) - datetime.timedelta(hours=self.timezone)
        year_difference = 2 - time_with_offset.year
        month_number = time_with_offset.month
        day = time_with_offset.day
        hour = time_with_offset.hour
        for year in self.years:
            declination = SPA(year - year_difference, month_number, day, hour, 0, 0, 0, self.latitude, self.longitude,
                              self.elevation, 0, 0, 0, 0)[0][3]
            declinations.append(declination)
        declination = sum(declinations) / len(declinations) * (1 if self.latitude >= 0 else -1)
        northwards_from_center = self.width / 2 * math.cos(math.radians(self.latitude)) * math.tan(math.radians(
            declination))
        return northwards_from_center

    def _foci(self) -> Tuple[float, float]:
        """
        Finds the foci of the ellipse.

        :return: The foci of the ellipse.

        c1, c2 = 0, 0
        F1 = - sqrt(a^2 - b^2) + c1, c2
        F2 = sqrt(a^2 - b^2) + c1, c2
        """
        foci_loc = math.sqrt((self.width / 2) ** 2 - (self.height / 2) ** 2)
        return -foci_loc, foci_loc

    def _angle(self, hour: NUMBER_TYPE) -> float:
        """
        Calculates the angle for a given hour.

        :param hour: The hour of the day.
        :return: The angle for a given hour.

        tan(θ) = tan(15° * hour) / sin(latitude)
        positive degrees = East, negative degrees = West
        0° = North, -90° = West, 90° = East
        θ = θ + angle_correction (if it is applied)
        """
        if self.correct_for_longitude:
            hour = hour + self._longitude_correction() / 60
        if hour <= 6:
            adj_hour = hour + 12
            flip = -1
        elif hour > 18:
            adj_hour = hour - 12
            flip = 1
        else:
            adj_hour = hour
            flip = 0
        hour_12 = adj_hour - 12
        for_tan = math.radians(15 * hour_12)
        tan = math.tan(for_tan)
        sin = math.sin(math.radians(self.latitude))
        tan_theta = tan / sin
        arc_tan = math.degrees(math.atan(tan_theta))
        arc = arc_tan
        if flip == -1:
            arc -= 180
        elif flip == 1:
            arc += 180
        if arc < -180:
            arc += 360
        elif arc > 180:
            arc -= 360
        return arc

    def _longitude_correction(self) -> NUMBER_TYPE:
        """
        Calculates the longitude correction.

        :return: The longitude correction in minutes.

        angle_correction = longitude - time_zone * 15
        minute_correction = angle_correction * 4
        """
        angle_correction = self.longitude - self.timezone * 15
        minute_correction = angle_correction * 4
        return minute_correction

    def _sundial(self, filename: Optional[str] = None) -> None:
        """
        Calculates data and draws the ellipse if filename is provided.

        :param filename: If filename is provided the sundial is saved.
        """
        ellipse = Ellipse((0, 0), self.width, self.height, fill=False)
        fig = plt.figure(1, dpi=200)
        ax = fig.add_subplot(1, 1, 1, aspect='equal')
        ax.fill(0, 0, alpha=0.2, facecolor='yellow',
                edgecolor='yellow', linewidth=1, zorder=1)

        ax.add_patch(ellipse)
        focis = self._foci()
        ax.plot([0], [self.height / 2], color='royalblue', marker='^', zorder=8)
        ax.plot([0], [-self.height / 2], color='royalblue', marker='v', zorder=7)
        ax.plot([self.width / 2], [0], color='royalblue', marker='>', zorder=6)
        ax.plot([-self.width / 2], [0], color='royalblue', marker='<', zorder=5)
        ax.plot([0], [0], color='k', marker='.', zorder=4)
        ax.scatter(list(focis), [0, 0], color='k', marker='.', zorder=3)
        ax.plot([-self.width / 2, self.width / 2], [0, 0], color='royalblue', zorder=1)
        ax.plot([0, 0], [-self.height / 2, self.height / 2], color='royalblue', zorder=2)

        zorder = -100

        def add_to_x(hour: NUMBER_TYPE) -> float:
            return (1 if self._angle(hour) > 0 else -1) * (1 if self.latitude > 0 else -1) - 0.1

        def add_to_y(hour: NUMBER_TYPE) -> float:
            return .7 if 90 > self._angle(hour) > -90 else -1

        for hour in range(24):
            arc = self._angle(hour)
            x_ans, y_ans = self._hour_location(arc)
            self.hour_locations[hour] = (arc, x_ans, y_ans)
            zorder -= 1
            ax.plot([x_ans], [y_ans], color='k', marker='.', zorder=zorder, linewidth=1)
            zorder -= 1
            ax.text(x_ans + add_to_x(hour) * 0.02 * self.width, y_ans + add_to_y(hour) * 0.02 * self.width, hour,
                    {'size': 4}, zorder=zorder)

        right = -1.
        zorder = -10
        for month in ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']:
            month_number = self.MONTH_TO_NUMBER[month]
            day = 1
            if month == 'Jul':
                right = .4
            zorder -= 1
            move_gnomon = self._move(month_number, day)
            self.gnomon_movement[f"{month} {day}"] = move_gnomon
            ax.scatter([0], [move_gnomon], color='mediumspringgreen', marker='_', zorder=zorder)
            zorder -= 1
            ax.text(right * self.width * 0.05, move_gnomon, f'{month} {day}', {'size': 4}, zorder=zorder)
            if month in ['Jun', 'Dec']:
                day = 21
                move_gnomon = self._move(month_number, day)
                self.gnomon_movement[f"{month} {day}"] = move_gnomon
                ax.scatter([0], [move_gnomon], color='mediumspringgreen', marker='_', zorder=zorder)
                ax.text(right * self.width * 0.05, move_gnomon, f'{month} {day}', {'size': 4}, zorder=zorder)

        if isinstance(filename, str):
            fig.savefig(filename)

    def _corrections(self, filename: Optional[str] = None) -> None:
        """
        Calculates the correction chart.

        :param filename: If filename is provided the correction chart is saved.
        """
        month_x = []
        month_y = []
        change_x = []
        change_y = []
        eots = []
        previous_value, previous_change = 0., 0.
        best_year = 1
        for year in self.years:
            month_points_x = []
            month_points_y = []
            change_points_x = []
            change_points_y = []
            eot_values = []
            is_leap = year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
            if best_year == 1 or not is_leap:
                best_year = year
            for days in range(1, 366 + int(is_leap)):
                day = datetime.datetime(year, 1, 1) + datetime.timedelta(days - 1)
                month = day.month
                day_of_month = day.day
                eot = SPA(year, month, day_of_month, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)[1][0]
                eot_values.append(eot + (self._longitude_correction() if not self.correct_for_longitude else 0))
            add_value = abs(min(eot_values)) + 10
            for index, eot_value in enumerate(eot_values):
                if index == 0:
                    previous_value = eot_values[0]
                    previous_change = -0.001
                    month_points_x.append(0)
                    month_points_y.append(previous_value)
                else:
                    change = ((eot_value + add_value) / (previous_value + add_value)) - 1
                    days = index + 1
                    day = datetime.datetime(year, 1, 1) + datetime.timedelta(days - 1)
                    month = day.month
                    day_of_month = day.day
                    if day_of_month == 1:
                        month_points_x.append(days - 1)
                        month_points_y.append(eot_value)
                    if change > 0 >= previous_change or change < 0 <= previous_change:
                        days = index  # index and not index + 1 because the max value happened the previous day
                        day = day - datetime.timedelta(1)
                        month = day.month
                        day_of_month = day.day
                        change_points_x.append(days - 1)
                        change_points_y.append(previous_value)
                    previous_change = change
                    previous_value = eot_value
            month_x.append(month_points_x)
            month_y.append(month_points_y)
            change_x.append(change_points_x)
            change_y.append(change_points_y)
            eots.append(eot_values)
        average_month_x = [sum(values) / len(values) for values in zip(*month_x)]
        average_month_y = [sum(values) / len(values) for values in zip(*month_y)]
        average_change_x = [sum(values) / len(values) for values in zip(*change_x)]
        average_change_y = [sum(values) / len(values) for values in zip(*change_y)]
        for index, value in enumerate(average_month_y):
            self.significant_eot[f"{self.NUMBER_TO_MONTH[index + 1]} 1"] = value
        for day_of_the_year, value in zip(average_change_x, average_change_y):
            date = datetime.datetime(best_year, 1, 1) + datetime.timedelta(day_of_the_year)
            self.significant_eot[f"{self.NUMBER_TO_MONTH[date.month]} {date.day}"] = value

        def order_function(item: Tuple[str, float]) -> int:
            return self.MONTH_TO_NUMBER[item[0].split()[0]] * 40 + int(item[0].split()[1])
        self.significant_eot = {k: v for k, v in sorted(self.significant_eot.items(), key=order_function)}
        self.average_eot = [sum(values) / len(values) for values in zip(*eots)]

        plt.figure(2)
        plt.plot(self.average_eot, zorder=0)
        plt.scatter(average_month_x, average_month_y, marker='.', color='r', zorder=1)
        plt.scatter(average_change_x, average_change_y, marker='.', color='g', zorder=2)

        if isinstance(filename, str):
            plt.savefig(filename)

    def create_sundial(self, sundial_filename: Optional[str] = None, corrections_filename: Optional[str] = None
                       ) -> None:
        """
        Creates the sundial and the correction chart.

        :param sundial_filename: If the filename is provided the sundial is saved.
        :param corrections_filename: If the filename is provided the correction chart is saved.
        """
        self._sundial(sundial_filename)
        self._corrections(corrections_filename)
        if self.show:
            plt.show()

    def __repr__(self) -> str:
        """
        Prints all the data.
        """
        str_repr = [f"Width: {round(self.width, 3)}, Height: {round(self.height, 3)}"]
        if self.correct_for_longitude:
            str_repr.append(f"Longitude correction: {round(self._longitude_correction(), 3)}")
        for date, value in self.gnomon_movement.items():
            str_repr.append(f"For {date}, the gnomon should move {round(value, 3)} meters forwards")
        for hour, values in self.hour_locations.items():
            str_repr.append(f"The angle for the {hour}th hour is {round(values[0], 3)} and the coordinates: {(round(values[1], 3), round(values[2], 3))}")
        for date, value in self.significant_eot.items():
            str_repr.append(f"For {date}, the equation of time (EOT) is {round(value, 3)} minutes")
        return "\n".join(str_repr)
