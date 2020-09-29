#!/usr/bin/env python

class BigForewordList:
    """ The list only consider the location after the current
    """
    def __init__(self, the_list):
        self.my_list = the_list
        self.pos = 0


    def filter(self, func):
        if self.my_list:
            return filter(func, self.my_list[self.pos:])
        else:
            return []


    def get_element(self, index):
        if not self.my_list:
            return None
        if index >= 0:
            if self.pos + index < len(self.my_list):
                return self.my_list[self.pos + index]
            else:
                return None
        else:
            if len(self.my_list) + index >= self.pos:
                return self.my_list[index]
            else:
                return None


    def len(self):
        return len(self.my_list) - self.pos


    def append(self, element):
        self.my_list.append(element)


    def get_current_element(self):
        if self.my_list and self.pos < len(self.my_list):
            return self.my_list[self.pos]
        else:
            return None


    def get_last_element(self):
        try:
            return self.my_list[-1]
        except IndexError:
            return None


    def move_foreword(self):
        if self.pos <= len(self.my_list) - 1:
            self.pos += 1
        if self.pos == 10000:
            del self.my_list[:10000]
            self.pos = 0


    def is_empty(self):
        return self.pos == len(self.my_list)


if __name__ == '__main__':
    print('Here be dragons...')