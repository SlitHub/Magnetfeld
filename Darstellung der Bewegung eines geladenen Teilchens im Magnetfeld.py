from math import sin, cos, tan, atan, asin, acos, pi, floor, atan2
from operator import add, sub
import pygame
import numpy as np

pygame.init()

# Constants, sizes and colours:

#Colour definitions as RGB code:

RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
PURPLE = (128, 0, 128)
ORANGE = (255, 165, 0)
GREY = (120, 120, 120)

# Defines the width and height of the window that's opened during execution:

DISPLAY_WIDTH = 1800
DISPLAY_HEIGHT = 1000

# Acts as conversion between meters and pixels. In this program one pixel equals 1/CORRECTOR meters. The "floor()" makes sure that only whole
# values can be displayed since the program crashes if no whole number (e.g. 1.5) is enetered. 

CORRECTOR = 0.0000005

CORRECTOR_X = floor(DISPLAY_WIDTH * 3 / 5)
CORRECTOR_Y = floor(DISPLAY_HEIGHT / 2)

POSITION_OF_FIELD = [DISPLAY_WIDTH*3/(CORRECTOR*5), DISPLAY_HEIGHT/(CORRECTOR*2), 0]


# Magnetic field constant and magnetic dipole moment of earth.

MAGNETIC_FIELD_CONSTANT = 4*pi*10**-7
MAGNETIC_DIPOLE_MOMENT = 6.496*10**22

# Mass and charge of relevant particles in Kilograms and Coulombs.

CHARGE_PROTON = 1.6021766208*10**-19
MASS_PROTON = 1.672621898*10**-29.9

CHARGE_ELECTRON = -1.6021766208*10**-19
MASS_ELECTRON = 9.10938356*10**-31

CHARGE_HELIUM = 2.0 * 1.6021766208*10**-19 
MASS_HELIUM = 4 * 1.660539040*10**-29.9


WINDOW = pygame.display.set_mode((DISPLAY_WIDTH, DISPLAY_HEIGHT))
pygame.display.set_caption("Magnetic Field Demo")
clock = pygame.time.Clock()


WINDOW.fill(WHITE)

"""
This function calculates the cross-product of two cartesian three component vectors.
"""

def crossproduct(u, v):
    return [u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]]
    


"""
These two functions allow me to freely transform between the cartesian and spherical/polar coordinate systems. This helps in terms of programming and imagination
because most forces will act perpendicular to the initial vector (cartesian system is handy) but the magnetic field is far easier to program in a spherical system.

As "cartesian" or "polar" a list with three components is expected. That is [x, y, z] for "cartesian" or [r, rho, phi]
(with both of the angles being in radiants) for "polar".

Basic trigonometry is used to achieve the transformation. The cartesian system takes 3 axis to pin one point in the coordinate system.
The spherical however takes a radius, one horizontal- and one equatorial angle.
"""
"""

def convert_to_polar(cartesian):
    polar = [0, 0, 0]
    # r = polar[0]
    polar[0] = ((cartesian[0])**2 + (cartesian[1])**2 + (cartesian[2])**2)**0.5
    # theta = polar[1]
    if cartesian[0] == 0 and cartesian[1] >= 0:
        polar[1] = pi/2
    elif cartesian[0] == 0 and cartesian[1] <= 0:
        polar[1] = 3/2*pi
    else:
        polar[1] = atan2(cartesian[0], cartesian[1])
    # phi = polar[2]
    if cartesian[2] == 0:
        polar[2] = pi/2
    else:
        polar[2] = atan2((cartesian[2]), (((cartesian[0])**2 + (cartesian[1])**2)**0.5))
    if polar[1] >= 2*pi:
        polar[1] = polar[1] % 2*pi 
    if polar[2] >= 2*pi:
        polar[2] = polar[2] % 2*pi 
    return polar

def convert_to_cartesian(polar):
    cartesian = [0, 0, 0]
    # x = cartesian[0]
    cartesian[0] = polar[0]*sin(polar[2])*cos(polar[1])
    # y = cartesian[1]
    cartesian[1] = polar[0]*sin(polar[2])*sin(polar[1])
    # z = cartesian[2]
    cartesian[2] = polar[0]*cos(polar[2])
    return cartesian
"""

def convert_to_polar(cartesian):
    (x, y, z) = cartesian
    r = (x**2 + y**2 + z**2)**0.5
    if r == 0 and y >= 0:
        theta = pi/2
    elif x == 0 and y <= 0:
        theta = 3/2*pi
    else:
        theta = atan2(x, y)
    if z == 0:
        phi = pi/2
    else:
        phi = atan2((x**2 + y**2)**0.5, z)
    if theta >= 2*pi:
        theta = theta % 2*pi
    if phi >= 2*pi:
        phi = phi % 2*pi
    return np.array((r, theta, phi))
        
def convert_to_cartesian(polar):
    (r, theta, phi) = polar
    x = r*sin(phi)*cos(theta)
    y = r*sin(phi)*sin(theta)
    z = r*cos(phi)
    return np.array((x, y, z))
        


"""
The center of the polar coordinate system is the center of the earth. This makes calculations far easier and helps with processing speed.
Because the particle will always be at a certain angle towards the center, said angle can also be easily extracted by simply converting the
cartesian coordinates to the polar ones with the function shown earlier. With the information of the particles position (which is predefined by the user),
the strength of the magnetic field at that position can be calculated. As for the magnetic field vector there are two points that need to be connected.
The length of the vector is the strength of the magnetic field that we already calculated. We are now can
calculate the force using the cross product (with the function that I already mentioned above) and we have the force vector. With the new force
vector and the mass of the respective particle we can calculate the acceleration in the force direction and with our time increments (time to the
next calculation) ultimately the new speed. One of the main problems this process is facing is that the particle changes direction at any
point in time thus there should be a new force calculated at the same rate.
"""

# Calculates the range from [0, 0, 0] to the magnetic field line with a given strength. The desired field strength and the azimuth angle is demandet.

def magnetic_field_radius(field_strength, azimuth_angle):
    return (((MAGNETIC_DIPOLE_MOMENT*MAGNETIC_FIELD_CONSTANT*(3*(sin(azimuth_angle))**2)**0.5)/(4*pi*field_strength))**(1/3))

# Calculates the strength of the magnetic field at the desired angle. A radius and the desired azimuth angle is expected.

def magnetic_field_strength(radius, azimuth_angle):
    return (MAGNETIC_DIPOLE_MOMENT*MAGNETIC_FIELD_CONSTANT*(3*(sin(azimuth_angle))**2)**0.5)/(4*pi*radius**3)

# This function draws a line that connects all points with the same magnetic field strength and the same equatorial angle. It takes a field_strength,
# the amount of increments (basically the resolution of the graph), the start x- and y-positions as pixels, the colour and the line thickness as arguments.
# It calculates the first point and adds a small number (2*pi/increments) to the azimuth angle. This will calculate a second point that serves as second
# statement for the pygame.draw.line() function. The function runs as long as the azimuth angle is not equal to 2*pi rad (360 deg).

def draw_magnetic_field(field_strength, increments, start_x, start_y, colour, line_thickness):
    azimuth_angle = 0
    while azimuth_angle <= (2*pi):
        coordinates_n0 = [sin(azimuth_angle)*magnetic_field_radius(field_strength, azimuth_angle), cos(azimuth_angle)*magnetic_field_radius(field_strength, azimuth_angle)]
        coordinates_n1 = [sin(azimuth_angle+2*pi/increments)*magnetic_field_radius(field_strength, azimuth_angle+2*pi/increments), cos(azimuth_angle+2*pi/increments)*magnetic_field_radius(field_strength, azimuth_angle+2*pi/increments)]
        azimuth_angle = azimuth_angle + (2*pi/increments)
        xy0 = ((coordinates_n0[0]*CORRECTOR + start_x), (coordinates_n0[1]*CORRECTOR + start_y))
        xy1 = ((coordinates_n1[0]*CORRECTOR + start_x), (coordinates_n1[1]*CORRECTOR + start_y))
        pygame.draw.line(WINDOW, colour, xy0, xy1, line_thickness)

# As with the draw_magnetic_field() function this one uses pretty much the sam principle: It converts the particle position to polar, takes its angles
# and calculates the field strength. Later it does the same but with a "0.0001" bigger azimuth angle which gives out two resulting polar vectors.
# Both of these are converted to cartesian and then subtracted from each other to leave a direction vector of the field. Each element of the vector is
# then divided by its sum and multiplied by the strength of the magnetic field at the particle location.
    
def magnetic_field_vector(particle_position_cart):
    (r, theta, phi) = convert_to_polar(particle_position_cart) #, (POSITION_OF_FIELD))
    b_strength = magnetic_field_strength(r, theta)
    b_radius_1 = magnetic_field_radius(b_strength, theta)
    b_radius_2 = magnetic_field_radius(b_strength, (theta + 0.001))
    direction_before = convert_to_cartesian([b_radius_1, theta, phi]) 
    direction_after = convert_to_cartesian([b_radius_2, (theta+0.001), phi])
    direction_raw = [direction_after[0] - direction_before[0], direction_after[1] - direction_before[1], direction_after[2] - direction_before[2]]
    return direction_raw / sum(direction_raw) * b_strength

# The speed vector is achieved by the input speed and the direction vector. Again, the vector is normed by dividing each element by the vectors sum and
# then multiplied by the flat speed. 

def particle_velocity(direction_vector, speed):
    return direction_vector / sum(direction_vector) * speed

# The lorentz force is essential for the particle movement. It calculates the force that is excerted on the particle. From Newton we know that F = m*a.
# It takes the crossproduct of the veloctiy and the magnetic field vector and multiplies each element of the list by the charge of the aprticle.

def lorentz_force(charge_particle, direction_vector, speed, particle_position):
    p_velocity = particle_velocity(direction_vector, speed)
    b_field_vector = magnetic_field_vector(particle_position)
    force_vector = np.cross(p_velocity, b_field_vector)
    return force_vector * charge_particle

# The particle can be made to move with adding/subtracting numbers of the list "particle_position". The list(map(add, x, y)) allows to easily and cleanly
# add all positions of two lists. The particle velocity is calculated, then the lorentz force and eventually the new position of the particle.
# Then a dot is drawn for each times the function runs. The count is simply to prevent an infinite loop. 

def particle_movement(charge_particle, direction_vector, speed, particle_position, particle_mass, colour):
    count = 10000
    for i in range(count):
        lorentz_vector = lorentz_force(charge_particle, direction_vector, speed, particle_position)
        acceleration = lorentz_vector / particle_mass
        particle_speed_normalized = particle_velocity(direction_vector,  speed)
        actual_speed = particle_speed_normalized + acceleration
        particle_position = particle_position + actual_speed
        particle_x = floor(particle_position[0]*CORRECTOR)
        particle_y = floor(particle_position[1]*CORRECTOR)
        pygame.draw.circle(WINDOW, colour, (particle_x, particle_y), 1, 0)
        yield particle_x, particle_y
        #particle_z = floor(particle_position[2]*CORRECTOR)
        

# LINES is a list of magnetic field lines that should be drawn. Components are the desired magnetic field strengths and the following loop runs
# "draw_magnetic_field()" for as many times as there are values in "LINES".

LINES = [0.00000025, 0.0000005, 0.00000075, 0.000001, 0.00000125, 0.0000015, 0.00000175, 2.469948941917352*10**-10.5, 2.469948941917352*10**-9, 2.469948941917352*10**-10, 2.469948941917352e-11]

for line in LINES:
    draw_magnetic_field(line, 360, floor(POSITION_OF_FIELD[0]*CORRECTOR), floor(POSITION_OF_FIELD[1]*CORRECTOR), BLUE, 1)

# This is the circle that depicts the earth. When the CORRECTOR and the
# display resolution are too small to draw the earth in its original size, the radius must be made higher (r must be greater or equal to 1).
# pygame.draw.circle(surface, colour, (x, y), radius, width)

#pygame.draw.circle(WINDOW, RED, (floor(DISPLAY_WIDTH/2), floor(DISPLAY_HEIGHT/2)), floor(6370000*CORRECTOR), 2)
pygame.draw.circle(WINDOW, RED, (floor(POSITION_OF_FIELD[0]*CORRECTOR), floor(POSITION_OF_FIELD[1]*CORRECTOR)), 1, 0)

# These are the various lines that get drawn once the program is running. As coordinates the inputs are display oriented. That mean with a Full HD
# display the particle will still start in the middle on the left. Currently the input in meters is (0, 71'428'571, 0)

particle_movement(CHARGE_ELECTRON, np.array([1, 0, 0]), 600000, np.array([0/CORRECTOR, ((1)*DISPLAY_HEIGHT)/(2*CORRECTOR), 0/CORRECTOR]), MASS_ELECTRON, GREEN)
particle_movement(CHARGE_PROTON, np.array([1, 0, 0]), 300000, np.array([0/CORRECTOR, ((1)*DISPLAY_HEIGHT)/(2*CORRECTOR), 0/CORRECTOR]), MASS_PROTON, PURPLE)
particle_movement(CHARGE_HELIUM, np.array([1, 0, 0]), 300000, np.array([0/CORRECTOR, ((1)*DISPLAY_HEIGHT)/(2*CORRECTOR), 0/CORRECTOR]), MASS_HELIUM, ORANGE)


running = True

print(magnetic_field_vector([20, 20, 20]))

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        
    pygame.display.update()
    clock.tick(30)

pygame.quit()
quit()


 
    