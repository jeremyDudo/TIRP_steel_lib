

# for plotting/saving
def subtitle(composition0):
    dependent_element = "Fe"
    plotName = dependent_element
    for key, value in sorted(composition0.items(), key=lambda item: item[1], reverse=True):
        plotName = plotName+"-%s%s" % (round(value*100,3), key)
    return plotName


# file name generator
def file_name(composition0, element1={}, element2={}):
    plotName = subtitle(composition0)
    
    if len(element1) == 0:
        return plotName

    title = "Vary " + element1["name"] + "-" + element2["name"] + "__(" + str(element1["start"]) + "-" + str(element1["end"]) + "-" + str(element1["length"]) + ")(" + str(element2["start"]) + "-" + str(element2["end"]) + "-" + str(element2["length"]) + ")"

    
    file_name = title + "_" + plotName 
    return file_name

