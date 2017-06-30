expected_str = "melon"
received_str = "apple"
basket = ["banana", "grapes", "strawberry", "melon", "orange"]
x = 0
step = int(raw_input("Input iteration step: "))
while(received_str != expected_str):
  if(x >= len(basket)): print "No more fruits left on the basket."; break
received_str = basket[x]
x += step # Change this to 3 to make the while statement
# evaluate to false, avoiding the break statement, using the else
clause
if(received_str==basket[2]): print "I hate",basket[2],"!"; break
if(received_str != expected_str): print "I am waiting for my
",expected_str,"."
else:
print "Finally got what I wanted! my precious ",expected_str,"!"
print "Going back home now !"
