DEFAULT_TEXT = "racecarenterelephantmalayalam"

def all_palindromes(text=DEFAULT_TEXT):
    """Return list with all palindrome strings within text.

    The base of this algorithm is to start with a given character,
    and then see if the surrounding characters are equal. If they
    are equal then it's a palindrome and is added to results set,
    extend to check if the next character on either side is equal,
    adding to set if equal, and breaking out of loop if not.

    This needs to be repeated twice, once for palindromes of
    odd lengths, and once for palindromes of an even length."""

    results = set()
    text_length = len(text)
    for idx, char in enumerate(text):

        # Check for longest odd palindrome(s)
        start, end = idx - 1, idx + 1
        while start >= 0 and end < text_length and text[start] == text[end]:
            results.add(text[start:end+1])
            start -= 1
            end += 1

        # Check for longest even palindrome(s)
        start, end = idx, idx + 1
        while start >= 0 and end < text_length and text[start] == text[end]:
            results.add(text[start:end+1])
            start -= 1
            end += 1

    return list(results)


if __name__ == '__main__':
    print all_palindromes()