import math
import numpy as np

"""

dotpitch = 0.282 # [mm/dot]
visual_angle = 10 # [degrees]
visual_range = 70 # [cm]

pixel = pixel_size(dotpitch,visual_angle,visual_range)
pixel = pixel * dotpitch 

dotpitch = 0.264 # [mm/dot]
visual_angle = 9.75 # [degrees]
visual_range = 60 # [cm]

pixel = pixel_size(dotpitch,visual_angle,visual_range)
pixel = pixel * dotpitch 

"""

def pixel_size(dotpitch,visual_angle,visual_range):
  
    # dotpitch = 0.282; #(mm) SMI
    # visual_angle = 0.29 #(deg)
    # visual_range = 60 #(cm)
    
    visual_angle = visual_angle * (math.pi / 180)
    a = visual_range * math.tan(visual_angle)
    a = a * 10
#    disp(concat([num2str(a),'mm']))
    pixel_num = a / dotpitch

    return pixel_num
    


"""
dotpitch = 0.264 # [mm/dot]
pixel_size = 80 / dotpitch
visual_distance = 60

pixel2angle(dotpitch,[pixel_size],visual_distance)


"""
def pixel2angle(dotpitch,pixel_num,visual_range):
    
    # if len(pixel_num) > 1:
    angle = []
    for pix in pixel_num:
        angle.append(math.atan(((pix*dotpitch)/10)/visual_range))
        
    ans = np.array(angle) * (180/math.pi)
    # else:
    # angle = math.atan(((pixel_num*dotpitch)/10)/visual_range)
    # ans = angle * (180/math.pi)
        
    return ans


def angle2cm(dotpitch,visual_angle,visual_range):
  
    # dotpitch = 0.282; #(mm) SMI
    # visual_angle = 0.29 #(deg)
    # visual_range = 60 #(cm)
    
    visual_angle = visual_angle * (math.pi / 180)
    a = visual_range * math.tan(visual_angle)
    a = a * 10
#    disp(concat([num2str(a),'mm']))
    pixel_num = a / dotpitch

    return pixel_num
    