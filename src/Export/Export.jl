# Please write a function that can export the entire EITFullData struct  into a file and save it at a specified location.
export save_eit_full_data, load_eit_full_data

using Serialization

function save_eit_full_data(data::EITFullData, filename::String)
    open(filename, "w") do io
        serialize(io, data)
    end
end


function load_eit_full_data(filename::String)
    open(filename, "r") do io
        return deserialize(io)
    end
end