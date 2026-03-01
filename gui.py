import matplotlib.pyplot as plt


class RealtimeTempPlot:
    """
    Рисует несколько кривых T(t) на одном графике.
    На вход обновления принимает температуру в Кельвинах, хранит/показывает в °C.
    """

    def __init__(self, title: str = "Температуры по участкам"):
        plt.ion()

        self.fig, self.ax = plt.subplots(figsize=(9, 5))
        self.ax.set_title(title)
        self.ax.set_xlabel("t, s")
        self.ax.set_ylabel("T, °C")
        self.ax.grid(True)

        self.series = {}
        self._legend_drawn = False

    def add_series(self, name: str, label: str | None = None, lw: float = 2.0):
        if label is None:
            label = name
        line, = self.ax.plot([], [], lw=lw, label=label)
        self.series[name] = {"t": [], "T": [], "line": line}
        self._legend_drawn = False

    def push_point(self, name: str, t: float, Tk: float):
        """Добавить точку. Tk — в Кельвинах, на графике будет в °C."""
        if name not in self.series:
            # авто-добавление серии, если забыли объявить
            self.add_series(name, label=name)

        Tc = Tk - 273.15
        s = self.series[name]
        s["t"].append(t)
        s["T"].append(Tc)
        s["line"].set_data(s["t"], s["T"])

    def redraw(self):
        if not self._legend_drawn:
            self.ax.legend()
            self._legend_drawn = True

        self.ax.relim()
        self.ax.autoscale_view()
        plt.pause(0.001)

    def close(self, block: bool = True):
        plt.ioff()
        plt.show(block=block)
        
    def hold(self):
        plt.ioff()
        plt.show()
        


def create_default_temp_plot():
    """
    Готовый набор линий под твой контур:
    выход АЗ, выход ПГ, выход опускного, конец контура
    """
    p = RealtimeTempPlot("T(t) в ключевых сечениях контура")
    p.add_series("AZ_in", label="АЗ (вход)")
    p.add_series("AZ_out", label="АЗ (выход)")
    p.add_series("PG_in", label="ПГ (вход)")
    p.add_series("PG_out", label="ПГ выход")
    return p

