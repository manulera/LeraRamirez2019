import matplotlib.pylab as plt
import numpy as np



def table_legend(handles,row_names, column_names,title_space=5,row_space=3,title_prespace=20):

    total_list =list()
    for l in handles:
        for ll in l:
            print "a"
            if type(ll)==list:
                total_list.append(ll[0])
            else:
                total_list.append(ll)



    title_space = " "*title_space

    title = title_prespace*" "+title_space.join(column_names)

    row_names = [c + row_space*" " for c in row_names]

    plt.legend(total_list, row_names + [''] * (len(total_list) - 2),
               title=title,
               ncol=len(column_names), numpoints=1, handletextpad=0, markerfirst=False)




x = np.linspace(0,1,20)
y1 = x
y2 = x*x

plt.figure()
a1=plt.plot(x,y1)
a2=plt.plot(x,y2)
b1=plt.scatter(x,y1)
b2=plt.scatter(x,y2)


table_legend([[a1,a2],[b1,b2]],["Line", "Parabola"],["Plot","Scatter"])
plt.show()