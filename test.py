from __future__ import division
from pylab import *
from extract import *

def test_find_objects():
    names = glob.glob("C:\\Users\\Nicholas\\Desktop\\1718+48\\*.fit")
    names.sort()
    photos = map(Photo, names)
    series = Series(photos)
    explore = Explore(search_rad = 10, pval=0.02, rad_fraction=0.75, extended_sum_rad=1)
    find = Find(explore=explore, check_dist=20)
    objects = series.find_objects(find=find)[0]
    f = 0.48764
    print names[int(len(photos) * f)]
    photo = photos[int(len(photos) * f)]
    photo.load()
    for obj in objects:
        if obj.states[0].time <= int(len(photos) * f) <= obj.states[-1].time:
            state = obj.track(int(len(photos) * f))
            position, radius = state.position, state.radius
            ((x_start, x_stop), (y_start, y_stop)) = box(position, radius, photo.scidata.shape)
            photo.scidata[x_start: x_stop, y_start: y_stop] *= 0.25
    photo.auto_plot()

def generate_lightcurve():
    names = glob.glob("C:\\Users\\Nicholas\\Desktop\\1718+48\\*.fit")
    names.sort()
    photos = map(Photo, names)
    series = Series(photos)
    explore = Explore(search_rad = 20, pval=0.02, rad_fraction=0.75, extended_sum_rad=1)
    find = Find(explore=explore, check_dist=20)
    found_objs = series.find_objects(find=find)[0]
    f = open("values.txt", "w")
    for i in range(0, len(found_objs), 10):
        objs = found_objs[i: i + 10]
        start_time = objs[0].states[0].time
        stop_time = objs[0].states[-1].time
        for time in range(start_time, stop_time + 1):
            photo = photos[time]
            print photo
            photo.load()
            line = str(photo) + " " + str(photo.headers["DATE"]) + ": "
            for obj in objs:
                state = obj.track(time)
                position, radius = state.position, state.radius
                intensity = photo.intensity(position, radius)
                line += str(intensity) + " "
            line = line[:-1] + "\n"
            f.write(line)
            photo.close()
    f.close()
    
def normalize_and_plot():
    f = open("values.txt")
    lines = f.readlines()
    f.close()
    times = []
    diffs = [[], [], [], [], [], [], [], [], []]
    for line in lines:
        split = line.split(" ")
        name_time = (split[1])[:-1]
        day_on = name_time.split("-")[-1]
        day = int(day_on.split("T")[0])
        hour_on = day_on.split("T")[1]
        hour = int(hour_on.split(":")[0])
        minute = int(hour_on.split(":")[1])
        second = float(hour_on.split(":")[2])
        time = second + 60*minute + 60*60*hour + 60*60*24*day
        times.append(time)
        for i in range(9):
            diff = diffs[i]
            diff.append(float(split[i + 2]) - float(split[i + 3]))
    avg = 15
    new_diffs = []
    for i in range(9):
        new_diffs.append([])
        for j in range(len(diffs[0]) - avg):
            new_diff = 0
            for s in range(j, j + avg):
                new_diff += diffs[i][s] / avg
            new_diffs[-1].append(new_diff)
    for i in range(3):
        plot(times[:-avg], array(new_diffs[i]) / 0.08)
    #show()
    for i in range(3, 6):
        plot(times[:-avg], array(new_diffs[i]) / 0.08)
    #show()
    for i in range(6, 9):
        plot(times[:-avg], array(new_diffs[i]) / 0.08)
    show()
    
if __name__=="__main__":
    normalize_and_plot()