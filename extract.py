from __future__ import division
from pylab import *
from scipy.stats import norm
from matplotlib import cm
from astropy.io import fits
from pybrain.tools.customxml.networkreader import NetworkReader
import shutil
import os
import sys
import ctypes
import platform
import glob

"""
Keyword Argument Classes:
"""

default_recommend_count = 10
default_recommend_maxiters = 100
default_recommend_max_rad = 25
default_recommend_min_rad = 3
default_recommend_rad_fraction = 0.8
default_recommend_space_fraction = 0.05
default_recommend_pval = 0.05
default_recommend_extended_sum_rad = 1

class Recommend:
    def __init__(self,
                 count=default_recommend_count,
                 maxiters=default_recommend_maxiters,
                 max_rad=default_recommend_max_rad,
                 min_rad=default_recommend_min_rad,
                 rad_fraction=default_recommend_rad_fraction,
                 space_fraction=default_recommend_space_fraction,
                 pval=default_recommend_pval,
                 extended_sum_rad=default_recommend_extended_sum_rad):
        """
        This class is a bundle of the keyword arguments for
        the recommend_stars function.
        """
        self.count = count
        self.maxiters = maxiters
        self.max_rad = max_rad
        self.min_rad = min_rad
        self.rad_fraction = rad_fraction
        self.space_fraction = space_fraction
        self.pval = pval
        self.extended_sum_rad = extended_sum_rad
    
    def __str__(self):
        """
        This method prints the keywords in
        an easy-to-read manner.
        """
        desc =  "Recommend object:\n"
        desc += "    count = " + str(self.count) + "\n"
        desc += "    maxiters = " + str(self.maxiters) + "\n"
        desc += "    max_rad = " + str(self.max_rad) + "\n"
        desc += "    min_rad = " + str(self.min_rad) + "\n"
        desc += "    rad_fraction = " + str(self.rad_fraction) + "\n"
        desc += "    space_fraction = " + str(self.space_fraction) + "\n"
        desc += "    pval = " + str(self.pval) + "\n"
        desc += "    extended_sum_rad = " + str(self.extended_sum_rad) + "\n"
        return desc
        
default_recommend = Recommend()

default_explore_search_rad = 30
default_explore_rad_fraction = 0.5
default_explore_pval = 0.05
default_explore_extended_sum_rad = 1

class Explore:
    def __init__(self,
                 search_rad=default_explore_search_rad,
                 rad_fraction=default_explore_rad_fraction,
                 pval=default_explore_pval,
                 extended_sum_rad=default_explore_extended_sum_rad):
        """
        This class is a bundle of the optional keyword arguments
        for the explore function.
        """
        self.search_rad = search_rad
        self.rad_fraction = rad_fraction
        self.pval = pval
        self.extended_sum_rad = extended_sum_rad
    
    def __str__(self):
        """
        This method prints the keywords in
        an easy-to-read manner.
        """
        desc =  "Explore object:\n"
        desc += "    search_rad = " + str(self.search_rad) + "\n"
        desc += "    rad_fraction = " + str(self.rad_fraction) + "\n"
        desc += "    pval = " + str(self.pval) + "\n"
        desc += "    extended_sum_rad = " + str(self.extended_sum_rad) + "\n"
        return desc
        
default_explore = Explore()

default_find_check_dist = 100
default_find_min_dim = 0.5
default_find_max_lap = 0.3
default_find_max_stretch = 20
default_find_min_breakable = 10

class Find:
    def __init__(self,
                 check_dist=default_find_check_dist,
                 min_dim=default_find_min_dim,
                 max_lap=default_find_max_lap,
                 max_stretch=default_find_max_stretch,
                 min_breakable=default_find_min_breakable,
                 recommend=default_recommend,
                 explore=default_explore):
        """
        This class is a bundle of the optional keyword arguments
        for the find_objects function.
        """
        self.check_dist = check_dist
        self.min_dim = min_dim
        self.max_lap = max_lap
        self.max_stretch = max_stretch
        self.min_breakable = min_breakable
        self.recommend = recommend
        self.explore = explore
    
    def __str__(self):
        """
        This method prints the keywords in
        an easy-to-use manner.
        """
        desc =  "Find object:\n"
        desc += "    check_dist = " + str(self.check_dist) + "\n"
        desc += "    min_dim = " + str(self.min_dim) + "\n"
        desc += "    max_lap = " + str(self.max_lap) + "\n"
        desc += "    max_stretch = " + str(self.max_stretch) + "\n"
        desc += "    min_breakable = " + str(self.min_breakable) + "\n"
        desc += "    recommend = Recommend object:\n"
        desc += "        count = " + str(self.recommend.count) + "\n"
        desc += "        maxiters = " + str(self.recommend.maxiters) + "\n"
        desc += "        max_rad = " + str(self.recommend.max_rad) + "\n"
        desc += "        min_rad = " + str(self.recommend.min_rad) + "\n"
        desc += "        rad_fraction = " + str(self.recommend.rad_fraction) + "\n"
        desc += "        space_fraction = " + str(self.recommend.space_fraction) + "\n"
        desc += "        pval = " + str(self.recommend.pval) + "\n"
        desc += "        extended_sum_rad = " + str(self.recommend.extended_sum_rad) + "\n"
        desc += "    explore = Explore object:\n"
        desc += "        search_rad = " + str(self.explore.search_rad) + "\n"
        desc += "        rad_fraction = " + str(self.explore.rad_fraction) + "\n"
        desc += "        pval = " + str(self.explore.pval) + "\n"
        desc += "        extended_sum_rad = " + str(self.explore.extended_sum_rad) + "\n"
        return desc

default_find = Find()

"""
Non-class Functions:
"""

def box((x, y), rad, shape, border=0):
    """
    This function returns:
        ((x_start, x_stop),
         (y_start, y_stop))
    given (x, y) and rad.
    If specified, it won't pick pixels 
    border away or closer to the walls.
    """
    x_start = max(x - rad, border)
    x_stop = min(x + rad, shape[0] - border)
    y_start = max(y - rad, border)
    y_stop = min(y + rad, shape[1] - border)
    return ((x_start, x_stop), (y_start, y_stop))

def get_free_space(dirname):
    """
    Returns folder/drive free space.
    I completely copied this function off the web.
    """
    if platform.system() == 'Windows':
        free_bytes = ctypes.c_ulonglong(0)
        ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(dirname),
                                                   None,
                                                   None,
                                                   ctypes.pointer(free_bytes))
        return free_bytes.value
    else:
        st = os.statvfs(dirname)
        return st.f_bavail * st.f_frsize

def copy(path_to_local, start, path_to_extern, available):
    """
    This function copies as much of the data as it can from "path_to_extern"
    to "path_to_local" starting at file number "start". The copied files won't
    occupy more than "available" bytes of the local available space.
    copy.status tells what fraction of the data has been copied.
    copy() returns a tuple of two items:
        The first item is whether all of the remaining data has been copied.
        The second item is, if not (1), what the next photo number to start on is.
        If (1), it is useless (it is len(photos))
    """
    copy.status = 0
    names = glob.glob(path_to_extern + "\\*.fit")
    names.sort()
    photo_size = os.path.getsize(names[start]) #all photos must have the same size.
    photos = min(available // photo_size, len(names) - start)
    for photonum in xrange(start, start + photos):
        src = names[photonum]
        dst = path_to_local + "\\" + src.split("\\")[-1]
        shutil.copyfile(src, dst)
        copy.status = (photonum - start + 1) / photos
        print copy.status
    return (available // photo_size >= len(names) - start, photonum + 1)
    
"""
Classes:
"""

class Photo:
    def __init__(self, name):
        """
        This function is the constructor of the Photo class.
        It takes the name of the file that will be loaded.
        """
        self.name = name

    def load(self):
        """
        This method loads the file.
        It creates self.hdulist and self.scidata
        It follows the (x, y) vs. (y, x) convention of MaxIm DL, not astropy.
        If there is an extra one layer thick dimension, that gets cut off.
        """
        self.hdulist = fits.open(self.name)
        self.scidata = self.hdulist[0].data
        self.scidata = self.scidata.transpose()
        if len(self.scidata.shape) == 3:
            blank_axis = list(self.scidata.shape).index(1)
            self.scidata = self.scidata.sum(axis=blank_axis)
        self.headers = self.hdulist[0].header

    def close(self):
        """
        This method deletes the contents of the opened file in RAM.
        """
        self.hdulist.close()
        try: del self.scidata
        except: pass
        try: del self.hdulist
        except: pass
        try: del self.headers
        except: pass
        """
        They may already be deleted,
        in which case we don't have to
        worry ourselves about them
        and can pass.
        """
    
    def __str__(self):
        """
        This function returns the name of the photo.
        """
        return self.name

    def intensity(self, (x, y), rad):
        """
        This method sums the pixels in the circle centered at (x, y) with radius rad.
        It also takes out backround noise according to the region around the circle.
        It approprietly corrects for stars near the edges of the photo.
        """
        inner_total = 0
        outer_total = 0
        inner_points = 0
        outer_points = 0
        ((x_start, x_stop), (y_start, y_stop)) = box((x, y), rad, self.scidata.shape)
        for xpos in xrange(x_start, x_stop):
            for ypos in xrange(y_start, y_stop):
                if (x - xpos)*(x - xpos) + (y - ypos)*(y - ypos) < rad*rad:
                    inner_total += self.scidata[xpos][ypos]
                    inner_points += 1
                else:
                    outer_total += self.scidata[xpos][ypos]
                    outer_points += 1
        intensity = max(inner_total - inner_points * outer_total / outer_points, 0)
        return intensity

    def stretch_plot(self, lower, upper):
        """
        stretch_plot plots the image.
        David helped me with the numpy minimally.
        """
        new_scidata = self.scidata.copy()
        new_scidata[new_scidata > upper] = upper
        new_scidata[new_scidata < lower] = lower
        new_scidata = (new_scidata - lower) / (upper - lower)
        imshow(new_scidata.transpose(), cmap = cm.Greys_r)

    def auto_plot(self):
        """
        auto_plot finds reasonable stretches for stretch_plot to perform.
        It then plots the image using that stretch.
        """
        data_pool = self.scidata.reshape(len(self.scidata)*len(self.scidata[0]))
        data_pool.sort()
        lower = data_pool[int(len(data_pool) * 0.00025)]
        upper = data_pool[int(len(data_pool) * 0.99975)]
        self.stretch_plot(lower, upper)

    def outer_stats(self, (x, y), rad):
        """
        outer_stats finds the mean and standard deviation of the pixels
        outside the circle of radius rad and center (x, y) but inside the box of
        side length 2*rad and center (x, y).
        This is taken as a measure of the background around a star of center
        (x, y) and radius rad.
        When at the border of the image, pixels are cut off approprietly.
        """
        ((x_start, x_stop), (y_start, y_stop)) = box((x, y), rad, self.scidata.shape)
        xvals, yvals = meshgrid(arange(x_start - x, x_stop - x),
                                arange(y_start - y, y_stop - y))
        xvals = xvals.transpose()
        yvals = yvals.transpose()
        mask = (xvals**2 + yvals**2 > rad**2)
        scidata_slice = (self.scidata[x_start: x_stop, y_start: y_stop])[mask]
        size = sum(mask)
        mean = sum(scidata_slice) / size
        std = sqrt(sum((scidata_slice - mean)**2) / (size - 1))
        return (mean, std)

    def recommend_stars(self, recommend=default_recommend):
        """
        This function recommends which stars are the best ones to pick.
        It is passed a Recommend instance which bundles all
        the optinal keyword arguments together.
        The first optional argumet, count, is the number of stars to find.
        Stars are picked as brighter by their brightest pixel.
        If a bright burst (not a star) is detected, it will not be selected.
        maxiters is the maximum number of searches possible
        (i.e. (maxiters - count) is the maximum number of bursts
        we'll tolerate before we just return what we've found)
        max_rad is the maximum radius a star can have.
        min_rad is the minimum radius a star can have.
        rad_fraction deals with how large we make the star's radius,
        while space_fraction deals with what fraction of the disk
        must be occupied to our bright pixels to be classified as a star.
        pval is the P-value cutoff of the criterion for being a star.
        extended_sum_rad is the number of pixels around the main pixel
        we average for our statistical analysis.
        Border stars are taken care of approprietly.
        David suggested I use slicing so long ago
        it isn't even apparent what part of the code
        I'm talking about.
        """
        stars = []
        allowed_mask = ones(self.scidata.shape)
        i = 0
        while len(stars) < recommend.count and i < recommend.maxiters:
            brightest_point = unravel_index((self.scidata * allowed_mask).argmax(), self.scidata.shape)
            ((x_start, x_stop), (y_start, y_stop)) = box(brightest_point,
                                                         recommend.max_rad,
                                                         self.scidata.shape,
                                                         border=recommend.extended_sum_rad + 1)
            outer_mean, outer_std = self.outer_stats(brightest_point,
                                                     recommend.max_rad + recommend.extended_sum_rad)
            mean_scidata_slice = []
            for x in xrange(x_start, x_stop):
                mean_scidata_slice.append([])
                for y in xrange(y_start, y_stop):
                    total_value = 0
                    size = 0
                    for x_sum in xrange(x - recommend.extended_sum_rad,
                                        x + recommend.extended_sum_rad + 1):
                        for y_sum in xrange(y - recommend.extended_sum_rad,
                                            y + recommend.extended_sum_rad + 1):
                            if allowed_mask[x_sum, y_sum]:
                                size += 1
                                total_value += self.scidata[x_sum, y_sum]
                    if size > 0:
                        mean_scidata_slice[-1].append(total_value / size)
                    else:
                        mean_scidata_slice[-1].append(0)
            allowed_slice = allowed_mask[x_start: x_stop, y_start: y_stop]
            cutoff = norm.ppf(1 - recommend.pval) / (2*recommend.extended_sum_rad + 1)
            #The divisor came from sqrt((2*recommend.extended_sum_rad + 1)**2).
            star_points = (mean_scidata_slice * allowed_slice > outer_mean + cutoff*outer_std)
            star_size = sum(star_points)
            if star_size > 0:
                xpos, ypos = meshgrid(arange(x_start, x_stop), arange(y_start, y_stop))
                xpos = xpos.transpose()
                ypos = ypos.transpose()
                center = (int(sum(xpos[star_points]) / star_size),
                          int(sum(ypos[star_points]) / star_size))
                rads = sqrt((xpos[star_points] - center[0])**2
                          + (ypos[star_points] - center[1])**2).astype(int)
                rads.sort()
                rad = rads[int(star_size * recommend.rad_fraction)] + 2
                fraction_occ = star_size / (pi*rad**2)
                if rad >= recommend.min_rad + 2 and fraction_occ > recommend.space_fraction:
                    stars.append((center, rad + 1))
                if len(stars) != recommend.count and i != recommend.maxiters - 1:
                    allowed_mask[xpos, ypos] = False
            else:
                return stars
            i += 1
        return stars

    def explore(self,
                center,
                allowed,
                explore=default_explore):
        """
        Given a star center of a previous photo,
        it finds the closest star center and radius in this photo.
        It must also be passed an argument of which pixels are valid to search over.
        The four optional arguments are:
            search_rad, or how large an area we search over for the star.
            rad_fraction, what fraction of the star-brightness disk we include in the star.
            pval, the P-value of our star-search.
            extended_sum_rad, the same as in recommend_stars().
        These are tied up into an instance of the Explore class.
        Stars at the border are taken care of approprietly.
        """
        ((x_start, x_stop), (y_start, y_stop)) = box(center, explore.search_rad, self.scidata.shape)
        old_scidata_slice = self.scidata[x_start: x_stop, y_start: y_stop]
        old_allowed_slice = allowed[x_start: x_stop, y_start: y_stop]
        brightest_point_relative = unravel_index((old_scidata_slice * old_allowed_slice).argmax(),
                                                 old_scidata_slice.shape)
        brightest_point = (brightest_point_relative[0] + x_start,
                           brightest_point_relative[1] + y_start)
        ((x_start, x_stop), (y_start, y_stop)) = box(brightest_point,
                                                     explore.search_rad,
                                                     self.scidata.shape,
                                                     border=explore.extended_sum_rad + 1)
        mean_scidata_slice = []
        for x in xrange(x_start, x_stop):
            mean_scidata_slice.append([])
            for y in xrange(y_start, y_stop):
                total_value = 0
                size = 0
                for x_sum in xrange(x - explore.extended_sum_rad,
                                    x + explore.extended_sum_rad + 1):
                    for y_sum in xrange(y - explore.extended_sum_rad,
                                        y + explore.extended_sum_rad + 1):
                        if allowed[x_sum, y_sum]:
                            size += 1
                            total_value += self.scidata[x_sum, y_sum]
                if size > 0:
                    mean_scidata_slice[-1].append(total_value / size)
                else:
                    mean_scidata_slice[-1].append(0)
        allowed_slice = allowed[x_start: x_stop, y_start: y_stop]
        outer_mean, outer_std = self.outer_stats(brightest_point, explore.search_rad)
        cutoff = norm.ppf(1 - explore.pval) / (2*explore.extended_sum_rad + 1)
        #The divisor came from sqrt((2*explore.extended_sum_rad + 1)**2).
        star_points = (mean_scidata_slice * allowed_slice > outer_mean + cutoff*outer_std)
        star_size = sum(star_points)
        if star_size > 0:
            xpos, ypos = meshgrid(arange(x_start - brightest_point[0],
                                         x_stop - brightest_point[0]),
                                  arange(y_start - brightest_point[1],
                                         y_stop - brightest_point[1]))
            xpos = xpos.transpose()
            ypos = ypos.transpose()
            new_center = (int(sum(xpos[star_points]) / star_size) + brightest_point[0],
                          int(sum(ypos[star_points]) / star_size) + brightest_point[1])
            rads = (sqrt((xpos[star_points])**2 + (ypos[star_points])**2)).astype(int)
            rads.sort()
            rad = rads[int(star_size * explore.rad_fraction)] + 2
            return (new_center, rad)
        else:
            print "There is no evidence of a star here any more!"
            return (center, int(explore.search_rad * explore.rad_fraction))
            
class State:
    def __init__(self, time, position, radius):
        """
        This class is comically simple.
        It is a placeholder for a star at a given:
            time (self.time)
            position (self.position)
            radius (self.radius)
        """
        self.time = time
        self.position = position
        self.radius = radius
    
    def __str__(self):
        """
        This function returns a string version
        of the State class.
        """
        message = "State(time = " + str(self.time) 
        message += ", position = " + str(self.position)
        message += ", radius = " + str(self.radius) + ")"
        return message
        
class Trackable:
    def __init__(self, states):
        """
        This class represents a star that can be tracked across time.
        It is passed a list of states the star is known to occupy.
        The states must be sorted by increasing time.
        """
        self.states = states

    def track(self, time):
        """
        This method locates the object at a given time.
        It returns an estimation of the object's state at time "time".
        time must be between the time of the first frame and the time
        of the last frame.
        """
        lower = 0
        upper = len(self.states) - 1
        while (upper - lower > 1):
            middle = (upper + lower) // 2
            if self.states[middle].time < time:
                lower = middle
            else:
                upper = middle
        start_time, stop_time = self.states[lower].time, self.states[upper].time
        start_pos, stop_pos = self.states[lower].position, self.states[upper].position
        start_rad, stop_rad = self.states[lower].radius, self.states[upper].radius
        try:
            time_position = (int((stop_pos[0]-start_pos[0]) * (time-start_time) / (stop_time-start_time) + start_pos[0]),
                             int((stop_pos[1]-start_pos[1]) * (time-start_time) / (stop_time-start_time) + start_pos[1]))
            time_radius = int((stop_rad-start_rad) * (time-start_time) / (stop_time-start_time) + start_rad)
            state = State(time, time_position, time_radius)
            return state
        except:
            state = State(start_time, start_pos, start_rad)
            return state
    
    def __str__(self):
        """
        This function returns a string
        version of the Trackable class.
        It consists of a list of the 
        string versions of the State classes.
        """
        message = "Trackable("
        for statenum, state in enumerate(self.states):
            message += str(state)
            if statenum != len(self.states) - 1:
                message += ", "
        return message + ")"

class Series:
    def __init__(self, photos):
        """
        This class represents a series of frames.
        It is what will be interfaced with from the outside
        (i.e. we will talk to the Series class when we want to detect exoplanets
        in the photographs taken)
        You pass in a list of the photographs (of class Photo)
        The photos shouldn't be loaded.
        """
        self.photos = photos
        
    def exoplanet_search(self,
                         find=default_find):
         """
         This method searches for exoplanets.
         The output will have the format:
             (exostar1_streak, exostar2_streak, ...)
         where an exostar is a star with an exoplanet, and a streak is
         a list of states in which the exostar was observed to have exoplanetary
         behaviour.
         At least 5 stars must be tracked.
         """
         stars, deleted = self.find_objects(find=find)
         print str(deleted / len(self.photos)) + "% of the data was ignored"
         """
         There must be an integer multiple of 5 stars
         in stars, and the stars must be grouped together in lumps
         of 5.
         """
         exostreaks = []
         net = NetworkReader.readFrom("../../Identifier/network.xml")
         for starnum in range(0, len(stars), 5):
             search_stars = stars[starnum: starnum + 5]
             start_time = search_stars[0].states[0].time
             stop_time = search_stars[0].states[-1].time
             for photonum in range(start_time, stop_time + 1, 10):
                 print self.photos[photonum]
                 photonum = min(photonum, stop_time - 10)
                 intensities = []
                 for slide in range(photonum, photonum + 10):
                     intensities.append([])
                     photo = self.photos[slide]
                     photo.load()
                     for star in search_stars:
                         state = star.track(slide)
                         brightness = photo.intensity(state.position, state.radius)
                         intensities[-1].append(brightness)
                     photo.close()
                 inpt = []
                 for starothernum in range(5):
                     lightcurve = []
                     for slides_from_zero in range(10):
                         lightcurve.append(intensities[slides_from_zero][starothernum])
                     array_version = array(lightcurve)
                     array_version /= average(array_version)
                     inpt += list(array_version)
                 nnet_output = net.activate(tuple(inpt))
                 for o in range(5):
                     if nnet_output[o] > 0.5:
                         exostreak = []
                         for slide in range(photonum, photonum + 10):
                             state = search_stars[o].track(slide)
                             exostreak.append(state)
                         exostreaks.append(exostreak)
         return exostreaks
 
    def find_objects_no_check(self,
                              photos,
                              start_time,
                              check_dist=100,
                              recommend=default_recommend,
                              explore=default_explore):
        """
        find_objects_no_check finds trackables
        which represent bright stars.
        These will later be searched for exoplanets.
        start_time is the time of the first photograph.
        check_freq is the distance between times we will calibrate states to
        encompass their stars.
        Recommend and Explore instances can also be passed as 
        optional keyword arguments - they will be propogated
        down the line to those functions when called if 
        specified. The tracker may fail - that
        is not this function's concern. 
        """
        i = 0
        stars = []
        while len(stars) < 10:
            photo = photos[i]
            print photo
            photo.load()
            stars = photo.recommend_stars(recommend=recommend)
            photo.close()
            i += 1
        i -= 1
        trackables = [Trackable([State(i + start_time, star[0], star[1])]) for star in stars]
        i += check_dist
        while i < len(photos):
            photo = photos[i]
            print photo
            photo.load()
            allowed = ones(photo.scidata.shape)
            for trackable in trackables: #sorted in decreasing order of brightness.
                last_state = trackable.states[-1]
                last_center = last_state.position
                new_center, new_radius = photo.explore(last_center,
                                                       allowed=allowed,
                                                       explore=explore)
                ((xmin, xmax), (ymin, ymax)) = box(new_center, new_radius, photo.scidata.shape)
                xpos, ypos = meshgrid(arange(xmin, xmax), arange(ymin, ymax))
                xpos.transpose()
                ypos.transpose()
                allowed[xpos, ypos] = False
                new_state = State(i + start_time, new_center, new_radius)
                trackable.states.append(new_state)
            photo.close()
            i += check_dist
        i -= check_dist
        if i != len(photos) - 1:
            i = len(photos) - 1
            photo = photos[i]
            print photo
            photo.load()
            allowed = ones(photo.scidata.shape)
            for trackable in trackables: #sorted in decreasing order of brightness.
                last_state = trackable.states[-1]
                last_center = last_state.position
                new_center, new_radius = photo.explore(last_center,
                                                       allowed=allowed,
                                                       explore=explore)
                ((xmin, xmax), (ymin, ymax)) = box(new_center, new_radius, photo.scidata.shape)
                xpos, ypos = meshgrid(arange(xmin, xmax), arange(ymin, ymax))
                xpos.transpose()
                ypos.transpose()
                allowed[xpos, ypos] = False
                new_state = State(i + start_time, new_center, new_radius)
                trackable.states.append(new_state)
            photo.close()
        return trackables                     
 
    def find_objects(self,
                     photos=None,
                     start_time=0,
                     find=default_find):
        """
        find_objects finds trackables which represent bright stars.
        It calls find_objects_no_check in binary-search fashion
        until it has a good list of trackables.
        It returns a tuple containing:
            (1) a list of trackables
            (2) the number of frames deleted
        """
        if photos is None:
            photos = self.photos
        print "len(photos) = " + str(len(photos))
        objects = self.find_objects_no_check(photos,
                                             start_time,
                                             find.check_dist,
                                             find.recommend,
                                             find.explore)
        if self.accept(objects, find.min_dim, find.max_lap, find.max_stretch):
            return (objects, 0)
        elif len(photos) <= find.min_breakable:
            print "Ignoring " + str(len(photos)) + " files."
            return ([], len(photos))
        else:
            prev_half = photos[:len(photos)//2]
            last_half = photos[len(photos)//2:]
            prev_objects, prev_deleted = self.find_objects(prev_half,
                                                           start_time,
                                                           find=find)
            last_objects, last_deleted = self.find_objects(last_half,
                                                           start_time + len(photos)//2,
                                                           find=find)
            return (prev_objects + last_objects, prev_deleted + last_deleted)
        
    def accept(self,
               objects,
               min_dim,
               max_lap,
               max_stretch):
        """
        This method determines whether 
        the objects returned by find_objects
        can be accepted.
        The three ways a trackable list 
        may not be accepted are:
            (1) The final intensities 
            are much lower than the initial
            intensities
            (2) two boxes are on top
            of each other.
            (3) the distances between
            the boxes has changed significantly.
        min_dim is the minimum fraction of the
        original star intensity still present. 
        max_lap is the maximum fraction of the area
        of one of two boxes that can overlap.
        max_stretch is the maximum distance change
        between boxes acceptible.
        """
        first_time = objects[0].states[0].time
        last_time = objects[0].states[-1].time
        first_states = [obj.track(first_time) for obj in objects]
        last_states = [obj.track(last_time) for obj in objects]
        first_photo = self.photos[first_time]
        first_photo.load()
        shape = first_photo.scidata.shape #This is used for calculating error 2
        first_intensities = array([first_photo.intensity(state.position,
                                                         state.radius) for state in first_states])
        first_photo.close()
        last_photo = self.photos[last_time]
        last_photo.load()
        last_intensities = array([last_photo.intensity(state.position,
                                                       state.radius) for state in last_states])
        last_photo.close()
        lost_mask = (last_intensities < first_intensities * min_dim)
        if any(lost_mask): print "error no. 1"; return False
        for state1num, state1 in enumerate(last_states):
            ((x1_start, x1_stop), (y1_start, y1_stop)) = box(state1.position,
                                                             state1.radius,
                                                             shape)
            for state2num, state2 in enumerate(last_states):
                if state1num != state2num:
                    ((x2_start, x2_stop), (y2_start, y2_stop)) = box(state2.position,
                                                                     state2.radius,
                                                                     shape)
                    if    x1_start <= x2_start <= x1_stop <= x2_stop: xlap = x1_stop - x2_start
                    elif  x1_start <= x2_start <= x2_stop <= x1_stop: xlap = x2_stop - x2_start
                    elif  x2_start <= x1_start <= x1_stop <= x2_stop: xlap = x1_stop - x1_start
                    elif  x2_start <= x1_start <= x2_stop <= x1_stop: xlap = x2_stop - x1_start
                    else: xlap = 0
                    if    y1_start <= y2_start <= y1_stop <= y2_stop: ylap = y1_stop - y2_start
                    elif  y1_start <= y2_start <= y2_stop <= y1_stop: ylap = y2_stop - y2_start
                    elif  y2_start <= y1_start <= y1_stop <= y2_stop: ylap = y1_stop - y1_start
                    elif  y2_start <= y1_start <= y2_stop <= y1_stop: ylap = y2_stop - y1_start
                    else: ylap = 0
                    lap = xlap*ylap
                    box_size = min((x1_stop - x1_start)*(y1_stop - y1_start),
                                   (x1_stop - x1_start)*(y2_stop - y2_start))
                    if lap >= box_size * max_lap: print "error no. 2"; return False
        first_distances = array([[sqrt((state1.position[0] - state2.position[0])**2
                                     + (state1.position[1] - state2.position[1])**2)
                                  for state1 in first_states]
                                 for state2 in first_states])
        last_distances = array([[sqrt((state1.position[0] - state2.position[0])**2
                                    + (state1.position[1] - state2.position[1])**2)
                                 for state1 in last_states]
                                for state2 in last_states])
        unchanged_mask = (abs(first_distances - last_distances) <= max_stretch)
        if not all(unchanged_mask): print "error no. 3"
        else: print "accepting"
        return all(unchanged_mask)