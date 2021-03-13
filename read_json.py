import json
import numpy as np
import matplotlib.pyplot as plt
    
############ JSON data layout #################
#there is a collection of objects called 'vstrands' which are the helices you place on the lattice; data are tracked as:
#       'row': row in the lattice
#       'col': col in the lattice
#       'num': number assigened to the scaffold layout panel
#       'scaf': a list of elements for the scaffold on the strand going from bp position 0 to #
#       'stap': a list of elements for the stpales on the strand going from bp position 0 to #
#       'loop': a track of which bp positions has a loop element; 0 for none, 1 for an insert
#       'skip': a track of which bp positions has a skip element; 0 for none, -1 for an insert
#       'scafLoop': ???
#       'stapLoop': ???
#       'stap_colors': a list of the colors for the staples along the loop, not sure how its assigned...  ???
#each element of 'scaf' or 'staple' is a list of the neighbors in the following order [helix_num_left, bp_position_left, helix_number_right, bp_position_right]
#if there is no staple of scaffold on a site the element is [-1,-1,-1,-1]. If you are at the end of a strand then there will be a [-1,-1, #,#] or [#,#,-1,-1]


def find_crossovers(data):
    #to find the crossover points we only need to search for elements in the 'stap' list that have neighbors on different strands
    #to know the scaffold element that the crossover occurs at we need to associated a position 'num', 'bp' to a scaffold element
    #for each element we can find how it connects to other ones and, once we have the starting position and end position can define a number that goes from 0 to
    #the length of the scaffold.

    vstrands = data['vstrands']

    crossovers = []

    num_tot = 0
    bp_long = len(vstrands[0]['stap'])

    #print(num_tot, bp_long)

    cnt_staple_bp = 0

    for vstrand in vstrands:
        for i in range(len(vstrand['stap'])):
            num = vstrand['num']
            if num+1>num_tot: num_tot = num+1 #keeps track of how many helices are in the structure

            stap = vstrand['stap'][i]
            if stap!=[-1,-1,-1,-1]:   #this occurs on elements where there is no staple strand
                cnt_staple_bp+=1

            if stap[0]!=stap[2]:  #there is a cross over if the neighboring elements are on different helices, the num is different
                if stap[0]>-1:    #at the end of a staple the element will be [-1, -1, X, X] or [X, X, -1, -1]
                    if stap[2]>-1:
                        crossovers.append([stap,[num,i]])

    print('num of strands:', num_tot)
    print('number of crossovers:', len(crossovers)/2) #crossovers are double counted

    #the crossover list is [[staple neighbors],[staple location]]
    #to find a staple crossover pair identify which strand the current crossover location is

    ###############################
    # find the locations of each end of the crossover staple
    staple_pairs = []

    for crossing in crossovers:
        next_tos = crossing[0]
        location = tuple(crossing[1])
        loc_left = (next_tos[0], next_tos[1])
        loc_right = (next_tos[2], next_tos[3])

        #two cases of which strand it could be
        if loc_left[0]==location[0]:
            paired_loc = loc_right
        elif loc_right[0]==location[0]:
            paired_loc = loc_left

        if (paired_loc,location) not in staple_pairs: #need to remove the duplications from this list
            staple_pairs.append((location, paired_loc))

    return staple_pairs, (num_tot, bp_long)
        
##########################################
      
def number_scaffold(data, params):
    #here we loop over the scaffold and number the order of elements, they are stored in a matrix
    #first we make a matrix of all the neighbor elements
    scaf_neighbors = np.empty(params, dtype=object)
    vstrands = data['vstrands']

    # we have different cases of the case where the scaffold has been broken or not
    flag_broken = False #checks to see if there is a break in the scaffold
    start_broken = np.s_[0,0]
    start_full = np.s_[0,0]

    neighbor_cnt = 0

    #to find the labeling of the scaffold first populate a martrix of neighbors while parsing json file
    for vstrand in vstrands:
        for i in range(len(vstrand['scaf'])):
            num = vstrand['num']
            scaf = vstrand['scaf'][i]

            if sum(scaf)>-4:
                scaf_neighbors[num,i] = scaf
                neighbor_cnt+=1
                start_full = np.s_[num, i]
                if scaf[0]+scaf[1]==-2:
                    flag_broken = True
                    start_broken = np.s_[num,i]

    print('scaffold count:', neighbor_cnt)
    #now we want to go through the array and keep track of the neighbors
    #while doing this fill a new array with the 
    scaf_position = np.zeros(params, dtype=int)

    if flag_broken: start_position = start_broken
    else: start_position = start_full

    #initalize the position
    scaf_position[start_position] = 1
    next_position = np.s_[scaf_neighbors[start_position][2], scaf_neighbors[start_position][3]]
    prev_position = start_position

    cnt = 0
    while (next_position!=start_position):
        current_position = next_position

        next_left, next_right = np.s_[scaf_neighbors[current_position][0], scaf_neighbors[current_position][1]], np.s_[scaf_neighbors[current_position][2], scaf_neighbors[current_position][3]]
        #print('next candidates:', next_left, next_right)

        if next_left == prev_position: next_position = next_right
        else: next_position = next_left

        if next_position == np.s_[-1,-1]: break

        #update the scaf_position array
        scaf_position[current_position] = scaf_position[prev_position]+1

        #print(prev_position, current_position, next_position)

        prev_position = current_position
        cnt+=1
    
    return scaf_position

####################################
    
#now we want to plot these crossovers in a nice plot
def distance_between(end_points, scaf_length):
    thetas = np.linspace(0,2*np.pi,scaf_length+1)
    th1, th2 = thetas[end_points[1]], thetas[end_points[0]]
    #shift one of these to th1 to zero
    v1 = (np.cos(th1), np.sin(th1))
    v2 = (np.cos(th2), np.sin(th2))
    distance = np.arccos(v1[0]*v2[0] + v2[1]*v1[1])/2./np.pi

    return distance

def plot_connection(end_points, x,y, axs, scaf_length):
    #given a pair of scaffold poistions (n1, n2) plot a line between them on a circle.
    #want the lines of the scaffold crossovers to be curved and have a curvation that is 0 when the distance between them is half the scaffold
    #and when the are close, when the curvature curvature to be half the distance between them.

    if end_points[0]==0: return -1
    if end_points[1]==0: return -1

    x1, x2 = x[end_points[0]-1], x[end_points[1]-1]
    y1, y2 = y[end_points[0]-1], y[end_points[1]-1]
    #axs.scatter([x1,x2], [y1,y2], c='k')

    def get_color(end_points, scaf_length):
        length = distance_between(end_points, scaf_length)
        R = 2*length
        return [R,0,0]

    axs.plot([x1,x2], [y1,y2], c=get_color(end_points, scaf_length))

def plot_circle_connection(end_points, x,y, axs, scaf_length):

    if end_points[0]==0: return -1
    if end_points[1]==0: return -1

    x1, x2 = x[end_points[0]-1], x[end_points[1]-1]
    y1, y2 = y[end_points[0]-1], y[end_points[1]-1]

    def get_color(end_points, scaf_length):
        length = distance_between(end_points, scaf_length)
        R = 2*length
        return [R,0,1-R]

    curvature = np.sqrt(4/((x1-x2)**2 + (y1-y2)**2)  - 1)
    if curvature==0:
        axs.plot([x1,x2], [y1,y2], c=get_color(end_points, scaf_length), lw=0.5)
        return

    else:
        L = 1/curvature

    s = np.sqrt((x1-x2)**2 + (y1-y2)**2)
    x3,y3 = (x1+x2)/2., (y1+y2)/2.

    xc_p = x3 + np.sqrt(L**2 - (s/2)**2)*(y1-y2)/s
    xc_m = x3 - np.sqrt(L**2 - (s/2)**2)*(y1-y2)/s

    yc_p = y3 + np.sqrt(L**2 - (s/2)**2)*(x2-x1)/s
    yc_m = y3 - np.sqrt(L**2 - (s/2)**2)*(x2-x1)/s

    #two possible centers, want the one that is outside the r=1 scaff circ
    d_p = xc_p**2 + yc_p**2
    d_m = xc_m**2 + yc_m**2

    #print(d_p, d_m)
    if d_p>d_m:
        center = (xc_p, yc_p)
    else:
        center = (xc_m, yc_m)

    #plot the circle
    th_0 = np.arctan2(center[1],center[0])
    th_0 = (th_0-np.pi)

    #opening angle
    v1, v2 = [x1-center[0], y1-center[1]], [x2-center[0], y2-center[1]]
    v_dot = v1[0]*v2[0] + v1[1]*v2[1]
    v_norm = np.sqrt(v1[0]*v1[0] + v1[1]*v1[1])*np.sqrt(v2[0]*v2[0] + v2[1]*v2[1])
    th_open = np.arccos(v_dot/v_norm)

    th_circ = np.linspace(th_0-th_open, th_0+th_open, 100)
    x_circ, y_circ = L*np.cos(th_circ)+center[0], L*np.sin(th_circ)+center[1]
    loc_inside = np.where(x_circ**2 + y_circ**2<1)

    axs.plot(x_circ[loc_inside], y_circ[loc_inside], c=get_color(end_points, scaf_length), lw=0.5)
    
def plot_connections(scaf_position, staple_pairs):
    
    #now for each staple pair we can find the scaffold positions that they connect
    scaf_connections = []
    for pair in staple_pairs:
        left_end, right_end = np.s_[pair[0][0], pair[0][1]], np.s_[pair[1][0], pair[1][1]]
        scaf_connections.append([scaf_position[left_end], scaf_position[right_end]])
        #print(scaf_connections[-1])

    thetas = np.linspace(0,2*np.pi,np.max(scaf_position)+1)
    xs = np.cos(thetas)
    ys = np.sin(thetas)

    fig, axs = plt.subplots(figsize=(8,8))

    plt.plot(np.cos(np.linspace(0,2*np.pi,100)),np.sin(np.linspace(0,2*np.pi,100)), 'k--')

    scaf_length = np.max(scaf_position)
    print('scaf_length', scaf_length)

    lengths = []

    cnt=0
    for connection in scaf_connections:
        plot_circle_connection(connection, xs, ys, axs, scaf_length)
        if connection[0]!=0:
            if connection[1]!=0:
                lengths.append(distance_between(connection,scaf_length))

    plt.figure(figsize=(4,3))
    plt.hist(np.array(lengths)*100., color='w', edgecolor='k')
    plt.xlim(0,50)
    #plt.ylim(0,250)
    plt.xlabel('scaffold distance [%]', fontsize=12)
    plt.ylabel('number of crossover', fontsize=12)
    plt.tight_layout()
    plt.show()


##############################################

if __name__ == "__main__":
    
    src = './16-helix-bundle-Ross/'
    filename = 'anaconda_v3.2.json'

    with open(src+filename,'r') as json_file:
        data = json.load(json_file)
    
    staple_pairs, lattice_params = find_crossovers(data)
    scaffold_position = number_scaffold(data, lattice_params)
    plot_connections(scaffold_position, staple_pairs)