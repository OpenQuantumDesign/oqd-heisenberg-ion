def standardize_string_format(string_input):

        temp = "".join([letter if letter.isalnum() else "_" for letter in string_input])

        if temp != string_input:
            return "".join([word.capitalize() for word in temp.strip().split("_")])
        else:
              return string_input