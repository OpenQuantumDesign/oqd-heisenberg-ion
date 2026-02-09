def convert_to_pascal(string_input):

    temp = "".join([letter if letter.isalnum() else "_" for letter in string_input])
    return "".join([word.capitalize() for word in temp.strip().split("_")])


def convert_to_snake_case(string_input):

    return "".join([letter.lower() if letter.isalnum() else "_" for letter in string_input])


def pascal_to_snake(string_input):

    result = []
    for i, c in enumerate(string_input):
        if c.isupper() and i > 0:
            if string_input[i - 1].islower():
                result.append("_")
            elif i + 1 < len(string_input) and string_input[i + 1].islower():
                result.append("_")

        result.append(c.lower())
    return "".join(result)
