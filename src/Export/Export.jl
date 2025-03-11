# Please write a function that can export the entire EITFullData struct  into a file and save it at a specified location.
export save_eit_full_data, load_eit_full_data


"""
Save the EITFullData struct to a file using Julia's serialization.

Arguments:
- data::EITFullData: The EITFullData object to be saved.
- filename::String: The path to the file where the data will be saved.

This function opens the file in write mode ("w"), serializes the EITFullData object, 
and writes it to the file in binary format.

Note: The file will be overwritten if it already exists.
"""
function save_eit_full_data(data::EITFullData, filename::String)
    open(filename, "w") do io
        serialize(io, data)
    end
end

"""
Load the EITFullData struct from a file using Julia's deserialization.

Arguments:
- filename::String: The path to the file from which to load the data.

Returns:
- EITFullData: The loaded EITFullData object.

This function opens the file in read mode ("r"), deserializes the data from the file, 
and returns the EITFullData object.

Note: The file must have been saved using save_eit_full_data or a similar serialization method.
"""
function load_eit_full_data(filename::String)
    open(filename, "r") do io
        return deserialize(io)
    end
end