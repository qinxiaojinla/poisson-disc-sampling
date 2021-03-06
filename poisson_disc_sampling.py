from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn
from math import sin,cos,sqrt,pi
import random
import numpy
import scipy.signal as signal
from itertools import product


class Grid:


    def __init__(self, r, *size):
        self.r = r

        self.size = size
        self.dim = len(size)

        self.cell_size = r/(sqrt(self.dim))

        self.widths = [int(size[k]/self.cell_size) + 1 for k in range(self.dim)]

        nums = product(*(range(self.widths[k]) for k in range(self.dim)))

        self.cells = {num:-1 for num in nums}
        self.samples = [] #已经采样的点
        self.active = [] #可选的采样点

    def reset(self):
        #初始化grid中的active points和sample points
        self.samples = []
        self.active = []

        for item in self.cells:
            self.cells[item] = -1

    def generate_new_points(self, point):
        #在(r,2r)之间的环面产生new points
        rad = random.uniform(self.r,2*self.r)
        angs = [random.uniform(0, 2*pi)]
        angs.extend(random.uniform(-pi/2, pi/2) for _ in range(self.dim-2))

        angs[0] = 2*angs[0]

        return self.convert(point, rad, angs)


    def poisson_disc(self, seed, k=30):
        #产生poisson-disc采样点
        self.reset()

        self.samples.append(seed)
        self.active.append(0)
        self.update(seed, 0)

        while self.active:

            idx = random.choice(self.active)
            point = self.samples[idx]
            new_point = self.genreate_points(k, point)

            if new_point:
                self.samples.append(tuple(new_point))
                self.active.append(len(self.samples)-1)
                self.update(new_point, len(self.samples)-1)
            else:
                self.active.remove(idx)

        return self.samples

    def genreate_points(self, k, point):
        #利用genreate产生K个new_point并利用这K个点用check函数检验point能不能成为sample point
        n = k
        while n:
            new_point = self.generate_new_points(point)
            if self.check(point, new_point):
                return new_point

            n -= 1

        return False

    def check(self, point, new_point):
        #计算point周围是否有neighbors，以及new_point的周围是否有neighbors
        #对于new_points的neighbors，计算距离是否大于r
        for i in range(self.dim):
            if not (0 < new_point[i] < self.size[i] or
                    self.point_in_which_cell(new_point) == -1):
                return False

        for item in self.find_neighbors(self.point_in_which_cell(point)):
            if self.point_distance(self.samples[item], new_point) < self.r**2:
                return False

        for item in self.find_neighbors(self.point_in_which_cell(new_point)):
            if self.point_distance(self.samples[item], new_point) < self.r**2:
                return False

        return True

    def convert(self, point, rad, angs):
        #把从以active point为坐标原点的极坐标系生成的随机点转换成直角坐标系下的点
        new_point = [point[0] + rad*cos(angs[0]), point[1] + rad*sin(angs[0])]
        if len(angs) > 1:
            new_point.extend(point[i+1] + rad*sin(angs[i]) for i in range(1,len(angs)))
        return new_point

    def point_in_which_cell(self, point):
        #返回point所在的cell
        return tuple(point[i]//self.cell_size for i in range(self.dim))

    def point_distance(self, tup1, tup2):
        #计算point之间的距离
        return sum((tup1[k] - tup2[k])**2 for k in range(self.dim))

    def cell_distance(self, tup1, tup2):
        #计算cell之间的距离是否小于2
        return sum(abs(tup1[k]-tup2[k]) for k in range(self.dim)) <= 2

    def find_neighbors(self, cell):
        #计算point的neighbors
        return (self.cells[tup] for tup in self.cells
                if self.cells[tup] != -1 and
                self.cell_distance(cell, tup))

    def update(self, point, index):
        #new_point产生后对grid进行更新
        self.cells[self.point_in_which_cell(point)] = index

    def __str__(self):
        return self.cells.__str__()

def data_processing(items):
    return ([item[i] for item in items] for i in range(len(items[0])))
