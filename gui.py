import matplotlib.pyplot as plt


class RealtimeTempPlot:

    def __init__(self, title="График", ylabel="Value"):
        plt.ion()

        self.fig, self.ax = plt.subplots(figsize=(9, 5))
        self.ax.set_title(title)
        self.ax.set_xlabel("t, часы")
        self.ax.set_ylabel(ylabel)
        self.ax.grid(True)

        self.series = {}
        self._legend_drawn = False

    def add_series(self, name: str, label=None, lw=2.0):
        if label is None:
            label = name

        line, = self.ax.plot([], [], lw=lw, label=label)

        self.series[name] = {
            "t": [],
            "Y": [],
            "line": line
        }

        self._legend_drawn = False

    def push_point(self, name: str, t: float, value: float):

        if name not in self.series:
            self.add_series(name)

        s = self.series[name]

        s["t"].append(t)
        s["Y"].append(value)

        s["line"].set_data(s["t"], s["Y"])

    def redraw(self):
        if not self._legend_drawn:
            self.ax.legend()
            self._legend_drawn = True

        self.ax.relim()
        self.ax.autoscale_view()
        plt.pause(0.001)

    def close(self, block=True):
        plt.ioff()
        plt.show(block=block)

    def hold(self):
        plt.ioff()
        plt.show()
        
def create_default_plot(parameter_name, unit, series):

    p = RealtimeTempPlot(
        title=f"{parameter_name}(t)",
        ylabel=f"{parameter_name}, {unit}"
    )

    for key, label in series.items():
        p.add_series(key, label=label)

    return p